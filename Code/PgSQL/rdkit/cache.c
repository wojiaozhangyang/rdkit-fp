// $Id$
//
//  Copyright (c) 2010, Novartis Institutes for BioMedical Research Inc.
//  All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met: 
//
//     * Redistributions of source code must retain the above copyright 
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following 
//       disclaimer in the documentation and/or other materials provided 
//       with the distribution.
//     * Neither the name of Novartis Institutes for BioMedical Research Inc. 
//       nor the names of its contributors may be used to endorse or promote 
//       products derived from this software without specific prior written
//       permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
#include <postgres.h>
#include <fmgr.h>
#include <utils/memutils.h>

#include "rdkit.h"
#include "cache.h"

#define MAGICKNUMBER  0xBEEC0DED
/*
 * Deconstructed values cache
 */

typedef enum EntryKind {
  MolKind,
  BfpKind,
  SfpKind,
  ReactionKind
} EntryKind;

typedef struct ValueCacheEntry {
  Datum toastedValue;
  EntryKind kind;
  
  union {
    struct {
      Mol *value;
      bytea *sign;
      CROMol mol;
    } mol;
    struct {
      Bfp *value;
      BfpSignature *sign;
      CBfp fp;
    } bfp;
    struct {
      Sfp *value;
      bytea *sign;
      CSfp fp;
    } sfp;
    struct {
      Reaction *value;
      bytea *sign;
      CChemicalReaction rxn;
    } reaction;
  } detoasted;
  
  struct ValueCacheEntry *prev;
  struct ValueCacheEntry *next;
} ValueCacheEntry;

#define NENTRIES (16)

typedef struct ValueCache {
  uint32 magickNumber;
  MemoryContext ctx;
  int32 nentries;
  ValueCacheEntry *head;
  ValueCacheEntry *tail;
  ValueCacheEntry *entries[NENTRIES];

  void (*resetOrig) (MemoryContext context);
  void (*deleteOrig) (MemoryContext context);
} ValueCache;

/*********** Managing LRU **********/
static void 
moveFirst(ValueCache *ac, ValueCacheEntry *entry)
{
  /*
   * delete entry form a list
   */
  Assert( entry != ac->head );
  
  if ( entry == ac->tail ) {
    Assert( entry->next == NULL );
    ac->tail = entry->prev;
    if ( ac->tail ) {
      ac->tail->next = NULL;
    }
    else {
        ac->head = NULL;
    }
  }
  else {
    entry->prev->next = entry->next;
    entry->next->prev = entry->prev;
  }

  /*
   * Install into head 
   */

  Assert( ac->head != NULL );
  Assert( ac->tail != NULL );

  entry->next = ac->head;
  entry->prev = NULL;
  ac->head->prev = entry;
  ac->head = entry;
}

#define DATUMSIZE(d)    VARSIZE_ANY(DatumGetPointer(d)) 

static int
cmpDatum(Datum a, Datum b)
{
  int32   la = DATUMSIZE(a);
  int32   lb = DATUMSIZE(b);

  if (la == lb) {
    return memcmp(DatumGetPointer(a), DatumGetPointer(b), la);
  }
  return (la > lb) ? 1 : -1;
}

static void
cleanupData(ValueCacheEntry *entry)
{
  pfree(DatumGetPointer(entry->toastedValue));
  
  switch(entry->kind) {
  case MolKind:
    if (entry->detoasted.mol.value) {
      pfree(entry->detoasted.mol.value);
    }
    if (entry->detoasted.mol.mol) {
      freeCROMol(entry->detoasted.mol.mol);
    }
    if (entry->detoasted.mol.sign) {
      pfree(entry->detoasted.mol.sign);
    }  
    break;
  case BfpKind:
    if (entry->detoasted.bfp.value) {
      pfree(entry->detoasted.bfp.value);
    }
    if (entry->detoasted.bfp.fp) {
      freeCBfp(entry->detoasted.bfp.fp);
    }
    if (entry->detoasted.bfp.sign) {
      pfree(entry->detoasted.bfp.sign);
    }  
    break;
  case SfpKind:
    if (entry->detoasted.sfp.value) {
      pfree(entry->detoasted.sfp.value);
    }
    if (entry->detoasted.sfp.fp) {
      freeCSfp(entry->detoasted.sfp.fp);
    }
    if (entry->detoasted.sfp.sign) {
      pfree(entry->detoasted.sfp.sign);
    }  
    break;
  case ReactionKind:
    if (entry->detoasted.reaction.value) {
      pfree(entry->detoasted.reaction.value);
    }
    if (entry->detoasted.reaction.rxn) {
      freeChemReaction(entry->detoasted.reaction.rxn);
    }
    if (entry->detoasted.reaction.sign) {
      pfree(entry->detoasted.reaction.sign);
    }  
    break;
  default:
    elog(ERROR, "Unknown kind: %d", entry->kind);
  }
  
  memset(entry, 0, offsetof(ValueCacheEntry, prev));
}

static void
makeEntry(ValueCache *ac, ValueCacheEntry *entry, Datum value, EntryKind kind)
{
  entry->toastedValue = (Datum)MemoryContextAlloc(ac->ctx, DATUMSIZE(value));
  entry->kind = kind;
  memcpy(DatumGetPointer(entry->toastedValue), DatumGetPointer(value),
	 DATUMSIZE(value));
}

static int
cmpEntry(const void *a, const void *b)
{
  return cmpDatum((*(ValueCacheEntry**)a)->toastedValue,
		  (*(ValueCacheEntry**)b)->toastedValue); 
}

/*********** Managing cache structure **********/

typedef struct CacheHolder {
  MemoryContext ctx;
  ValueCache *cache;
  struct CacheHolder *next;
} CacheHolder;

static CacheHolder *holder = NULL;

static void
cleanupRDKitCache(MemoryContext context)
{
  CacheHolder *h = holder, *p = NULL;

  /*
   * Find holder and clean non-postgres values.
   * Note, one context could contains several caches
   */
  while (h) {
    if (h->ctx == context) {
      if (h->cache->ctx != context || h->cache->magickNumber != MAGICKNUMBER) {
	elog(WARNING, "Something wrong in cleanupRDKitCache");
      }
      else {
	int i;
	
	for (i=0;i<h->cache->nentries;i++) {
	  cleanupData(h->cache->entries[i]);
	}
	h->cache->nentries = 0;
      }

      /* remove current holder from list */
      if (p==NULL) {
	holder = h->next;
	free(h);
	h = holder;
      }
      else {
	p->next = h->next;
	free(h);
	h = p->next;
      }
      continue;
    }
    
    p = h;
    h = h->next;
  }
}

MemoryContextMethods *methodsOrig = NULL;
MemoryContextMethods methodsCache;

static void
resetCacheContext(MemoryContext context)
{
  cleanupRDKitCache(context);
  methodsOrig->reset(context);
}

static void
deleteCacheContext(MemoryContext context)
{
  cleanupRDKitCache(context);
  methodsOrig->delete_context(context);
}

static ValueCache*
createCache(void *cache, struct MemoryContextData * ctx) 
{
  ValueCache *ac;
  CacheHolder *newholder;

  if (cache != NULL) {
    ac = (ValueCache*)cache;
    if (ac->ctx != ctx) {
      elog(ERROR, "We can't use our approach with cache :(");
    }
  }

  /*
   * We need to add cleanup data to delete and reset of out memory context,
   * for that we install new handlers, but we need to store old ones.
   * HACK!: we believe that there is only single memory context type 
   * in postgres 
   */

  /* define new methods */
  if (!methodsOrig) {
    methodsOrig = ctx->methods;
    methodsCache = *methodsOrig;
    methodsCache.reset = resetCacheContext;
    methodsCache.delete_context = deleteCacheContext;
  }

  /*
   * Try to connect to existing cache!
   */
  newholder = holder;
  while(newholder) {
    if (newholder->ctx == ctx) {
      cache = (void*)newholder->cache;
      break;
    }
    newholder = newholder->next;
  }

  if (!cache) {
    /*
     * did not found a cache in current context, make new one
     */
    cache = MemoryContextAllocZero(ctx, sizeof(ValueCache));
    ac = (ValueCache*) cache;
    ac->magickNumber = MAGICKNUMBER;
    ac->ctx = ctx;
    
    newholder = malloc(sizeof(*newholder));
    if (!newholder) {
      elog(ERROR, "Could not allocate %ld bytes", sizeof(*newholder)); 
    }
    
    /* init holder */
    newholder->ctx = ctx;
    newholder->cache = ac;
    
    if (!(ctx->methods == methodsOrig ||
	  ctx->methods == &methodsCache /* already used for another cache */)) {
        elog(ERROR, "We can't use our approache with cache :((");
    }
    
    ctx->methods = &methodsCache;

    /* store holder */
    newholder->next = holder;
    holder = newholder;
  }

  return (ValueCache*)cache;
}

/*********** SEARCHING **********/

/*
 * entry->kind == MolKind
 *      m, mol, fp, val
 * entry->kind == FpKind
 *      f, fp, val
 */

static void
fetchData(ValueCache *ac, ValueCacheEntry *entry,
          void **detoasted, void **internal, void **sign)
{
  MemoryContext   old;
  void * _tmp;
        
  switch(entry->kind) {
  case MolKind:
    if (detoasted) {
      if (entry->detoasted.mol.value == NULL) {
	Mol *detoastedMol;
	
	detoastedMol = DatumGetMolP(entry->toastedValue);
	entry->detoasted.mol.value
	  = MemoryContextAlloc( ac->ctx, VARSIZE(detoastedMol));
	memcpy(entry->detoasted.mol.value, detoastedMol, VARSIZE(detoastedMol));
      }
      *detoasted = entry->detoasted.mol.value;
    }

    if (internal) {
      if (entry->detoasted.mol.mol == NULL) {
	fetchData(ac, entry, &_tmp, NULL, NULL);
	entry->detoasted.mol.mol = constructROMol(entry->detoasted.mol.value);
      }
      *internal = entry->detoasted.mol.mol;
    }

    if (sign) {
      if (entry->detoasted.mol.sign == NULL) {
	fetchData(ac, entry, NULL, &_tmp, NULL);
	old = MemoryContextSwitchTo( ac->ctx );
	entry->detoasted.mol.sign = makeMolSignature(entry->detoasted.mol.mol);
	MemoryContextSwitchTo(old);
      }
      *sign = entry->detoasted.mol.sign;
    }
    break;
  case BfpKind:
    if (detoasted) {
      if (entry->detoasted.bfp.value == NULL) {
	Bfp *detoastedFP;
	
	detoastedFP = DatumGetBfpP(entry->toastedValue);
	entry->detoasted.bfp.value
	  = MemoryContextAlloc( ac->ctx, VARSIZE(detoastedFP));
	memcpy(entry->detoasted.bfp.value, detoastedFP, VARSIZE(detoastedFP));
      }
      *detoasted = entry->detoasted.bfp.value;
    }

    if (internal) {
      if (entry->detoasted.bfp.fp == NULL) {
	fetchData(ac, entry, &_tmp, NULL, NULL);
	entry->detoasted.bfp.fp
	  = constructCBfp(entry->detoasted.bfp.value);
      }
      *internal = entry->detoasted.bfp.fp;
    }

    if (sign) {
      if (entry->detoasted.bfp.sign == NULL) {
	fetchData(ac, entry, NULL, &_tmp, NULL);
	old = MemoryContextSwitchTo( ac->ctx );
	entry->detoasted.bfp.sign
	  = makeBfpSignature(entry->detoasted.bfp.fp);
	MemoryContextSwitchTo(old);
      }
      *sign = entry->detoasted.bfp.sign;
    }
    break;
  case SfpKind:
    if (detoasted) { 
      if (entry->detoasted.sfp.value == NULL) {
	Sfp *detoastedFP;
	
	detoastedFP = DatumGetSfpP(entry->toastedValue);
	entry->detoasted.sfp.value
	  = MemoryContextAlloc( ac->ctx, VARSIZE(detoastedFP));
	memcpy(entry->detoasted.sfp.value, detoastedFP, VARSIZE(detoastedFP));
      }
      *detoasted = entry->detoasted.sfp.value;
    }

    if (internal) {
      if (entry->detoasted.sfp.fp == NULL) {
	fetchData(ac, entry, &_tmp, NULL, NULL);
	entry->detoasted.sfp.fp
	  = constructCSfp(entry->detoasted.sfp.value);
      }
      *internal = entry->detoasted.sfp.fp;
    }

    if (sign) {
      if (entry->detoasted.sfp.sign == NULL) {
	fetchData(ac, entry, NULL, &_tmp, NULL);
	old = MemoryContextSwitchTo( ac->ctx );
	entry->detoasted.sfp.sign
	  = makeSfpSignature(entry->detoasted.sfp.fp, NUMBITS);
	MemoryContextSwitchTo(old);
      }
      *sign = entry->detoasted.sfp.sign;
    }
    break;
  case ReactionKind:
    if (detoasted) {
      if (entry->detoasted.reaction.value == NULL) {
	Reaction *detoastedRxn;
	
	detoastedRxn = DatumGetReactionP(entry->toastedValue);
	entry->detoasted.reaction.value
	  = MemoryContextAlloc( ac->ctx, VARSIZE(detoastedRxn));
	memcpy(entry->detoasted.reaction.value, detoastedRxn,
	       VARSIZE(detoastedRxn));
      }
      *detoasted = entry->detoasted.reaction.value;
    }
    
    if (internal) {
      if (entry->detoasted.reaction.rxn == NULL) {
	fetchData(ac, entry, &_tmp, NULL, NULL);
	entry->detoasted.reaction.rxn
	  = constructChemReact(entry->detoasted.reaction.value);
      }
      *internal = entry->detoasted.reaction.rxn;
    }
    
    if (sign) {
      if (entry->detoasted.reaction.sign == NULL) {
	fetchData(ac, entry, NULL, &_tmp, NULL);
	old = MemoryContextSwitchTo( ac->ctx );
	entry->detoasted.reaction.sign
	  = makeReactionSign(entry->detoasted.reaction.rxn);
	MemoryContextSwitchTo(old);
      }
      *sign = entry->detoasted.reaction.sign;
    }
    break;
  default:
    elog(ERROR, "Unknown kind: %d", entry->kind);
  }
}

static void*
SearchValueCache(void *cache, struct MemoryContextData * ctx, 
		 /*  in: */ Datum a, EntryKind kind, 
		 /* out: */ void **detoasted, void **internal, void **sign)
{
  ValueCache *ac = createCache(cache, ctx);
  ValueCacheEntry *entry;

  /*
   * Fast check of recent used value 
   */
  if (ac->head && cmpDatum(ac->head->toastedValue, a) == 0) {
    Assert( ac->head->kind == kind );
    fetchData(ac, ac->head, detoasted, internal, sign);
    return cache;
  }

  if (ac->head == NULL) {
    ac->entries[0]
      = ac->head
      = ac->tail
      = MemoryContextAllocZero(ctx, sizeof(ValueCacheEntry));
    ac->nentries = 1;
    makeEntry(ac, ac->head, a, kind);
    fetchData(ac, ac->head, detoasted, internal, sign);
    return cache;
  }

  do {
    ValueCacheEntry **StopLow = ac->entries;
    ValueCacheEntry **StopHigh = ac->entries + ac->nentries;
    ValueCacheEntry **StopMiddle;
    int cmp;
    
    while (StopLow < StopHigh) {
      StopMiddle = StopLow + ((StopHigh - StopLow) >> 1);
      entry = *StopMiddle;
      cmp = cmpDatum(entry->toastedValue, a); 
      
      if (cmp == 0) {
	moveFirst(ac, entry);
	Assert( ac->head->kind == kind );
	fetchData(ac, ac->head, detoasted, internal, sign);
	return cache;
      }
      else if (cmp < 0) {
        StopLow = StopMiddle + 1;
      }
      else {
        StopHigh = StopMiddle;
      }
    }
  } while(0);

  /*
   * Not found 
   */

  if (ac->nentries < NENTRIES) {
    entry
      = ac->entries[ac->nentries]
      = MemoryContextAllocZero(ctx, sizeof(ValueCacheEntry));

    /* install first */
    entry->next = ac->head;
    entry->prev = NULL;
    ac->head->prev = entry;
    ac->head = entry;
    
    ac->nentries ++;
    
    makeEntry(ac, ac->head, a, kind);
    fetchData(ac, ac->head, detoasted, internal, sign);
  } 
  else {
    cleanupData(ac->tail);
    moveFirst(ac, ac->tail);
    makeEntry(ac, ac->head, a, kind);
    fetchData(ac, ac->head, detoasted, internal, sign);
  }

  qsort(ac->entries, ac->nentries, sizeof(ValueCacheEntry*), cmpEntry);   
  return cache;
}

void* 
searchMolCache(void *cache, struct MemoryContextData * ctx, Datum a, 
               Mol **m, CROMol *mol, bytea ** val)
{
  return  SearchValueCache(cache, ctx, 
                           /*  in: */ a, MolKind, 
                           /* out: */ (void**)m, (void**)mol, (void**)val);
}

void* 
searchBfpCache(void *cache, struct MemoryContextData * ctx, Datum a, 
	       Bfp **f, CBfp *fp, BfpSignature **val)
{
  return  SearchValueCache(cache, ctx, 
                           /*  in: */ a, BfpKind, 
                           /* out: */ (void**)f, (void**)fp, (void**)val);
}

void* 
searchSfpCache(void *cache, struct MemoryContextData * ctx, Datum a, 
	       Sfp **f, CSfp *fp, bytea **val)
{
  return  SearchValueCache(cache, ctx, 
                           /*  in: */ a, SfpKind, 
                           /* out: */ (void**)f, (void**)fp, (void**)val);
}

void* 
searchReactionCache(void *cache, struct MemoryContextData * ctx, Datum a, 
		    Reaction **rxn, CChemicalReaction *crxn, bytea ** val)
{
  return  SearchValueCache(cache, ctx, 
                           /*  in: */ a, ReactionKind, 
                           /* out: */ (void**)rxn, (void**)crxn, (void**)val);
}
