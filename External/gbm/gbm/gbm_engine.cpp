//  GBM by Greg Ridgeway  Copyright (C) 2003

#include "stdafx.h"
#include "gbm_engine.h"

CGBM::CGBM()
{
    adFadj = NULL;
    adZ = NULL;
    afInBag = NULL;
    aiNodeAssign = NULL;
    aNodeSearch = NULL;
    
    cDepth = 0;
    cMinObsInNode = 0;
    dBagFraction = 0.0;
    dLambda = 0.0;
    fInitialized = false;
    cTotalInBag = 0;
    cTrain = 0;
    cValid = 0;

    pData = NULL;
    pDist = NULL;
    pNodeFactory = NULL;
    ptreeTemp = NULL;
}


CGBM::~CGBM()
{
    if(adFadj != NULL)
    {
        delete [] adFadj;
        adFadj = NULL;
    }
    if(adZ != NULL)
    {
        delete [] adZ;
        adZ = NULL;
    }
    if(afInBag != NULL)
    {
        delete [] afInBag;
        afInBag = NULL;
    }
    if(aiNodeAssign != NULL)
    {
        delete [] aiNodeAssign;
        aiNodeAssign = NULL;
    }
    if(aNodeSearch != NULL)
    {
        delete [] aNodeSearch;
        aNodeSearch = NULL;
    }
    if(ptreeTemp != NULL)
    {
        delete ptreeTemp;
        ptreeTemp = NULL;
    }
    // must delete the node factory last!!! at least after deleting trees
    if(pNodeFactory != NULL)
    {
        delete pNodeFactory;
        pNodeFactory = NULL;
    }
}


GBMRESULT CGBM::Initialize
(
    CDataset *pData,
    CDistribution *pDist,
    double dLambda,
    unsigned long cTrain,
    double dBagFraction,
    unsigned long cDepth,
    unsigned long cMinObsInNode,
	bool takeAllPos
)
{
    GBMRESULT hr = GBM_OK;
    unsigned long i=0;

    if(pData == NULL)
    {
        hr = GBM_INVALIDARG;
        goto Error;
    }
    if(pDist == NULL)
    {
        hr = GBM_INVALIDARG;
        goto Error;
    }

    this->pData = pData;
    this->pDist = pDist;
    this->dLambda = dLambda;
    this->cTrain = cTrain;
    this->dBagFraction = dBagFraction;
    this->cDepth = cDepth;
    this->cMinObsInNode = cMinObsInNode;
	this->takeAllPos = takeAllPos ;

    // allocate the tree structure
    ptreeTemp = new CCARTTree;
    if(ptreeTemp == NULL)
    {
        hr = GBM_OUTOFMEMORY;
        goto Error;
    }

    cValid = pData->cRows - cTrain;
	   
	cTotalInBag = (unsigned long)(dBagFraction*cTrain);

    adZ = new double[cTrain];
    if(adZ == NULL)
    {
        hr = GBM_OUTOFMEMORY;
        goto Error;
    }
    adFadj = new double[pData->cRows];
    if(adFadj == NULL)
    {
        hr = GBM_OUTOFMEMORY;
        goto Error;
    }

    pNodeFactory = new CNodeFactory();
    if(pNodeFactory == NULL)
    {
        hr = GBM_OUTOFMEMORY;
        goto Error;
    }
    hr = pNodeFactory->Initialize(cDepth);
    if(GBM_FAILED(hr))
    {
        goto Error;
    }
    ptreeTemp->Initialize(pNodeFactory);

    // array for flagging those observations in the bag
    afInBag = new bool[cTrain];
    if(afInBag==NULL)
    {
        hr = GBM_OUTOFMEMORY;
        goto Error;
    }
    // aiNodeAssign tracks to which node each training obs belongs
    aiNodeAssign = new ULONG[cTrain];
    if(aiNodeAssign==NULL)
    {
        hr = GBM_OUTOFMEMORY;
        goto Error;
    }
    // NodeSearch objects help decide which nodes to split
    aNodeSearch = new CNodeSearch[2*cDepth+1];
    if(aNodeSearch==NULL)
    {
        hr = GBM_OUTOFMEMORY;
        goto Error;
    }
    for(i=0; i<2*cDepth+1; i++)
    {
        aNodeSearch[i].Initialize(cMinObsInNode);
    }
    vecpTermNodes.resize(2*cDepth+1,NULL);

    fInitialized = true;

Cleanup:
    return hr;
Error:
    goto Cleanup;
}




GBMRESULT CGBM::Predict
(
    unsigned long iVar,
    unsigned long cTrees,
    double *adF,
    double *adX,
    unsigned long cLength
)
{
    GBMRESULT hr = GBM_OK;


    return hr;
}


GBMRESULT CGBM::Predict
(
    double *adX,
    unsigned long cRow,
    unsigned long cCol,
    unsigned long cTrees,
    double *adF
)
{
    GBMRESULT hr = GBM_OK;


    return hr;
}



GBMRESULT CGBM::GetVarRelativeInfluence
(
    double *adRelInf,
    unsigned long cTrees
)
{
    GBMRESULT hr = GBM_OK;
    int iVar=0;

    for(iVar=0; iVar<pData->cCols; iVar++)
    {
        adRelInf[iVar] = 0.0;
    }

    return hr;
}


GBMRESULT CGBM::PrintTree()
{
    GBMRESULT hr = GBM_OK;

    hr = ptreeTemp->Print();
    if(GBM_FAILED(hr)) goto Error;

Cleanup:
    return hr;
Error:
    goto Cleanup;
}




GBMRESULT CGBM::iterate
(
    double *adF,
    double &dTrainError,
    double &dValidError,
    double &dOOBagImprove,
    int &cNodes
)
{
    GBMRESULT hr = GBM_OK;
    unsigned long i = 0;
    unsigned long cBagged = 0;
    int npos = 0;

    if(!fInitialized)
    {
        hr = GBM_FAIL;
        goto Error;
    }

    dTrainError = 0.0;
    dValidError = 0.0;
    dOOBagImprove = 0.0;

    vecpTermNodes.assign(2*cDepth+1,NULL);

    // randomly assign observations to the Bag
    cBagged = 0;
    npos = 0 ;
    for(i=0; i<cTrain; i++)
    {
 
		if((this->takeAllPos && this->pData->adY[i] == 1.0) || unif_rand()*(cTrain-i) < cTotalInBag-cBagged)
        {
            afInBag[i] = true;
            cBagged++;
			if (this->pData->adY[i] == 1.0)
				npos ++ ;
        }
        else
        {
            afInBag[i] = false;
        }
    }

    #ifdef NOISY_DEBUG
    fprintf(stderr,"Compute working response\n");
    #endif
    hr = pDist->ComputeWorkingResponse(pData->adY, 
                                       pData->adMisc,
                                       pData->adOffset,
                                       adF, 
                                       adZ,
                                       pData->adWeight,
                                       afInBag,
                                       cTrain);
    if(GBM_FAILED(hr))
    {
        goto Error;
    }

    #ifdef NOISY_DEBUG
    fprintf(stderr,"Reset tree\n");
    #endif
    hr = ptreeTemp->Reset();
    #ifdef NOISY_DEBUG
    fprintf(stderr,"grow tree\n");
    #endif
    hr = ptreeTemp->grow(adZ,pData,pData->adWeight,adFadj,
                         cTrain,cTotalInBag,dLambda,cDepth,
                         cMinObsInNode,
                         afInBag,
                         aiNodeAssign,aNodeSearch,vecpTermNodes);
    if(GBM_FAILED(hr))
    {
        goto Error;
    }

    #ifdef NOISY_DEBUG
    fprintf(stderr,"get node count\n");
    #endif
    hr = ptreeTemp->GetNodeCount(cNodes);
    if(GBM_FAILED(hr))
    {
        goto Error;
    }

    // Now I have adF, adZ, and vecpTermNodes (new node assignments)
    // Fit the best constant within each terminal node
    #ifdef NOISY_DEBUG
    fprintf(stderr,"fit best constant\n");
    #endif
    hr = pDist->FitBestConstant(pData->adY,
                                pData->adMisc,
                                pData->adOffset,
                                pData->adWeight,
                                adF,
                                adZ,
                                aiNodeAssign,
                                cTrain,
                                vecpTermNodes,
                                (2*cNodes+1)/3, // number of terminal nodes
                                cMinObsInNode,
                                afInBag,
                                adFadj);
    if(GBM_FAILED(hr))
    {
        goto Error;
    }

    // update training predictions
    // fill in missing nodes where N < cMinObsInNode
    hr = ptreeTemp->Adjust(aiNodeAssign,adFadj,cTrain,
                           vecpTermNodes,cMinObsInNode);
    if(GBM_FAILED(hr))
    {
        goto Error;
    }
    ptreeTemp->SetShrinkage(dLambda);

    dOOBagImprove = pDist->BagImprovement(pData->adY,
                                          pData->adMisc,
                                          pData->adOffset,
                                          pData->adWeight,
                                          adF,
                                          adFadj,
                                          afInBag,
                                          dLambda,
                                          cTrain);

    // update the training predictions
    for(i=0; i < cTrain; i++)
    {
        adF[i] += dLambda*adFadj[i];
    }
    dTrainError = pDist->Deviance(pData->adY,
                                  pData->adMisc,
                                  pData->adOffset,
                                  pData->adWeight,
                                  adF,
                                  cTrain);
    // update the validation predictions
    hr = ptreeTemp->PredictValid(pData,cValid,adFadj);
    for(i=cTrain; i < cTrain+cValid; i++)
    {
        adF[i] += adFadj[i];
    }
    if(pData->fHasOffset)
    {
        dValidError = 
            pDist->Deviance(&(pData->adY[cTrain]),
                            &(pData->adMisc[cTrain]),
                            &(pData->adOffset[cTrain]),
                            &(pData->adWeight[cTrain]),
                            &(adF[cTrain]),
                            cValid);
    }
    else
    {
        dValidError = pDist->Deviance(&(pData->adY[cTrain]),
                                      &(pData->adMisc[cTrain]),
                                      NULL,
                                      &(pData->adWeight[cTrain]),
                                      &(adF[cTrain]),
                                      cValid);
    }

Cleanup:
    return hr;
Error:
    goto Cleanup;
}


GBMRESULT CGBM::TransferTreeToRList
(
    int *aiSplitVar,
    double *adSplitPoint,
    int *aiLeftNode,
    int *aiRightNode,
    int *aiMissingNode,
    double *adErrorReduction,
    double *adWeight,
    double *adPred,    
    VEC_VEC_CATEGORIES &vecSplitCodes,
    int cCatSplitsOld
)
{
    GBMRESULT hr = GBM_OK;

    hr = ptreeTemp->TransferTreeToRList(pData,
                                        aiSplitVar,
                                        adSplitPoint,
                                        aiLeftNode,
                                        aiRightNode,
                                        aiMissingNode,
                                        adErrorReduction,
                                        adWeight,
                                        adPred,
                                        vecSplitCodes,
                                        cCatSplitsOld,
                                        dLambda);

    return hr;
}


