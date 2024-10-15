#include <cstdint>
#include <iostream>
#include <cmath>
#include <cstring>
#include <cstdlib>
#include <random>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <ctime>
#include <iomanip>
#include <vector>
#include "MurmurHash3.h"
#include <fstream>
#include <string>

#define HASH_SEED 51318471
//#define HASH_SEED 74151457
using namespace std;

string avgF1RatioFileName = "./resultsFreeRS/F1OfRatios.txt";
ofstream avgF1RatioFile;

struct SSKeyNode;

struct SSValueNode {
    uint32_t counter;
    SSValueNode *next;
    SSValueNode *pre;
    SSKeyNode *firstKeyNode;
};

struct SSKeyNode {
    uint32_t key;
    uint32_t err;
    SSValueNode *parent;
    SSKeyNode *next;
};

class StreamSummary {
private:
    uint32_t SS_MAX_SIZE;
    uint32_t size = 0;
    SSValueNode *firstValueNode = nullptr;
    unordered_map<uint32_t, SSKeyNode *> hashTable;
public:
    StreamSummary(uint32_t SS_MAX_SIZE) {
        this->SS_MAX_SIZE=SS_MAX_SIZE;
    }
    ~StreamSummary(){
        SSValueNode *ssValueNode = firstValueNode;
        while (ssValueNode != nullptr) {
            SSKeyNode *ssKeyNode = ssValueNode->firstKeyNode;
            do {
                SSKeyNode *tempKeyNode = ssKeyNode;
                ssKeyNode = ssKeyNode->next;
                free(tempKeyNode);
            } while (ssKeyNode != ssValueNode->firstKeyNode);

            SSValueNode *tempValueNode = ssValueNode;
            ssValueNode = ssValueNode->next;
            free(tempValueNode);
        }
    }

    uint32_t getSize() const{
        return size;
    }
    uint32_t getMinVal() const{
        if (size != 0) {
            return firstValueNode->counter;
        } else {
            return -1;
        }
    }

    SSValueNode* getFirstValueNode() const{
        return firstValueNode;
    }

    SSKeyNode* getMinKeyNode() const{
        return firstValueNode->firstKeyNode;
    }
    SSKeyNode* findKey(uint32_t key) const{
        auto iter = hashTable.find(key);
        if (iter == hashTable.end()) {
            return nullptr;
        } else {
            return iter->second;
        }
    }

    void SSPush(uint32_t key, uint32_t value,uint32_t err) {
        SSKeyNode *newKeyNode;
        SSValueNode *newValueNode= nullptr;
        if(size>=SS_MAX_SIZE){//when stream summary is full, drop a key node of the first value node (the smallest value node)
            //delete the second key node (or the first key node when the first value node only has one key node)
            SSKeyNode *ssKeyNode = firstValueNode->firstKeyNode->next;
            auto iter = hashTable.find(ssKeyNode->key);
            hashTable.erase(iter);

            if(ssKeyNode->next==ssKeyNode){//the first value node only has one key node
                newValueNode=firstValueNode;
                firstValueNode= firstValueNode->next;//may be nullptr
                firstValueNode->pre= nullptr;
            }else{
                firstValueNode->firstKeyNode->next = ssKeyNode->next;
            }
            size--;
            newKeyNode=ssKeyNode;
            newKeyNode->err=err;
        }else{
            newKeyNode = (SSKeyNode *) calloc(1, sizeof(SSKeyNode));
            newKeyNode->err=err;
        }

        SSValueNode *lastValueNode= nullptr;
        SSValueNode *curValueNode=firstValueNode;
        while(curValueNode!= nullptr and curValueNode->counter<value){//find the first value node larger than value
            lastValueNode=curValueNode;
            curValueNode=curValueNode->next;
        }

        if(curValueNode== nullptr or curValueNode->counter>value)//need to create a new value Node between lastValueNode and curValueNode
        {
            if(newValueNode== nullptr){
                newValueNode = (SSValueNode *) calloc(1, sizeof(SSValueNode));
                newValueNode->firstKeyNode = newKeyNode;
            }

            newKeyNode->key=key;
            newKeyNode->parent = newValueNode;
            newKeyNode->next = newKeyNode;

            newValueNode->counter = value;
            newValueNode->pre = lastValueNode;
            newValueNode->next=curValueNode;

            if(curValueNode!= nullptr){
                curValueNode->pre=newValueNode;
            }

            if(lastValueNode== nullptr){//means newValueNode should be the firstValueNode
                firstValueNode = newValueNode;
            }else{
                lastValueNode->next=newValueNode;
            }

            hashTable[key] = newKeyNode;
            size++;
        }else{//value==curValueNode->counter
            //add into curValueNode
            //the key nodes form a circular list
            //compared to inserting the new key node as the first key node, inserting the new key node after the first key node can avoid check whether the ori first key node's next is itself
            newKeyNode->key=key;
            newKeyNode->parent = curValueNode;
            newKeyNode->next = curValueNode->firstKeyNode->next;
            curValueNode->firstKeyNode->next = newKeyNode;

            hashTable[key] = newKeyNode;
            size++;
        }
    }
    //newVal must be larger than the old val
    void SSUpdate(SSKeyNode *ssKeyNode, uint32_t newValue) {

        SSValueNode *newValueNode= nullptr;
        SSValueNode *lastValueNode= nullptr;
        SSValueNode *curValueNode= nullptr;

        if (ssKeyNode->next == ssKeyNode) {//the value node only has one key node
            newValueNode=ssKeyNode->parent;
            lastValueNode=newValueNode->pre;
            curValueNode=newValueNode->next;

            if(lastValueNode!= nullptr){
                lastValueNode->next=curValueNode;
            }
            if(curValueNode!= nullptr){
                curValueNode->pre=lastValueNode;
            }
            if(firstValueNode==newValueNode){
                firstValueNode=newValueNode->next;//may be nullptr
                firstValueNode->pre= nullptr;
            }
        }
        else{
            lastValueNode=ssKeyNode->parent;
            curValueNode=lastValueNode->next;
            //detach the key node from old value node
            //Since the linked list is single linked, we have to traverse the list to find the pre node.
            //We choose to replace the node content instead of finding the pre node.
            SSKeyNode* nextKeyNode=ssKeyNode->next;
            uint32_t key=ssKeyNode->key;
            uint32_t err=ssKeyNode->err;
            ssKeyNode->key=nextKeyNode->key;
            nextKeyNode->key=key;
            ssKeyNode->err=nextKeyNode->err;
            nextKeyNode->err=err;

            hashTable[ssKeyNode->key]=ssKeyNode;
            hashTable[nextKeyNode->key]=nextKeyNode;

            ssKeyNode->next=nextKeyNode->next;
            if(nextKeyNode==lastValueNode->firstKeyNode){
                lastValueNode->firstKeyNode=ssKeyNode;
            }
            ssKeyNode=nextKeyNode;
        }

        while(curValueNode!= nullptr and curValueNode->counter<newValue){//find the first value node larger than newValue
            lastValueNode=curValueNode;
            curValueNode=curValueNode->next;
        }

        if(curValueNode== nullptr or curValueNode->counter>newValue)//need to create a new value Node between lastValueNode and curValueNode
        {
            if(newValueNode== nullptr){
                newValueNode = (SSValueNode *) calloc(1, sizeof(SSValueNode));
                newValueNode->firstKeyNode = ssKeyNode;
                ssKeyNode->parent = newValueNode;
                ssKeyNode->next = ssKeyNode;
            }

            newValueNode->counter = newValue;
            newValueNode->pre = lastValueNode;
            newValueNode->next=curValueNode;

            if(curValueNode!= nullptr){
                curValueNode->pre=newValueNode;
            }

            if(lastValueNode== nullptr){//means newValueNode should be the firstValueNode
                firstValueNode = newValueNode;
            }else{
                lastValueNode->next=newValueNode;
            }

        }else{//newValue==curValueNode->counter
            ssKeyNode->parent = curValueNode;
            ssKeyNode->next = curValueNode->firstKeyNode->next;
            curValueNode->firstKeyNode->next = ssKeyNode;

            if(newValueNode!= nullptr){
                free(newValueNode);
            }
        }
    }
    void getKeyVals(unordered_map<uint32_t,uint32_t>& keyVals) const{
        for(auto iter=hashTable.begin();iter!=hashTable.end();iter++){
            uint32_t key=iter->first;
            SSKeyNode* keyNode=iter->second;
            uint32_t value=keyNode->parent->counter;
            keyVals[key]=value;
        }
    }
};

class SKETCH{
private:
    uint32_t SS_MAX_SIZE;
    uint32_t REG_NUM;
    StreamSummary ss;
    uint8_t* regArray;
    double updateProb=1.0;//regs被一个新元素更新的概率

    inline uint32_t getRegValue(uint32_t regIdx){
        uint32_t regStartByteIdx=regIdx*5/8;
        uint32_t regStartBitIdx=regIdx*5%8;
        return ((*(uint16_t*)(regArray+regStartByteIdx))>>regStartBitIdx)&0x1F;
    }
    inline void setRegValue(uint32_t regIdx,uint32_t newValue){
        uint32_t regStartByteIdx=regIdx*5/8;
        uint32_t regStartBitIdx=regIdx*5%8;

        uint16_t* value_ptr=(uint16_t*)(regArray+regStartByteIdx);
        *value_ptr=(*value_ptr)&(~(((uint16_t)0x1F)<<regStartBitIdx));
        *value_ptr=(*value_ptr)|(newValue<<regStartBitIdx);
    }

public:
    SKETCH(uint32_t SS_MAX_SIZE,uint32_t REG_NUM):ss(SS_MAX_SIZE),REG_NUM(REG_NUM){
        this->SS_MAX_SIZE=SS_MAX_SIZE;
        regArray=new uint8_t[int(ceil(REG_NUM*5.0/8.0))]();
    }
    SKETCH(const SKETCH& c):SKETCH(c.SS_MAX_SIZE,c.REG_NUM){
        uint32_t byteNum=int(ceil(REG_NUM*5.0/8.0));
        for(uint32_t byteIdx=0;byteIdx<byteNum;byteIdx++){
            regArray[byteIdx]=c.regArray[byteIdx];
        }
    }
    void insert(const pair<uint32_t,uint32_t>& pkt){
        uint32_t hashValue1=0;
        MurmurHash3_x86_32(&pkt,8,HASH_SEED, &hashValue1);

        uint32_t hashValue2=0;
        MurmurHash3_x86_32(&pkt, 8,HASH_SEED+15417891, &hashValue2);

        uint32_t regIdx=hashValue1%REG_NUM;
        uint32_t rhoValue=min(31,__builtin_clz(hashValue2)+1);

        uint32_t oriRhoValue= getRegValue(regIdx);

        if(rhoValue>oriRhoValue){
            double updateValue=1.0/updateProb;
            if(rhoValue!=31){
                updateProb-=(pow(2,-1.0*oriRhoValue)-pow(2,-1.0*rhoValue))/REG_NUM;
            }else{//when a reg's value is 31, it cannot be updated again
                updateProb-=pow(2,-1.0*oriRhoValue)/REG_NUM;
            }
            setRegValue(regIdx,rhoValue);

            uint32_t updateValueToInt= floor(updateValue);

            double temp=(double)rand()/RAND_MAX;
            if(temp<(updateValue- updateValueToInt)){
                updateValueToInt+=1;
            }
            SSKeyNode* keyNode=ss.findKey(pkt.first);
            if(keyNode!= nullptr){
                uint32_t newVal=updateValueToInt+keyNode->parent->counter;
                ss.SSUpdate(keyNode,newVal);
                // cout<<"newVal"<<newVal<<endl;
            }else{
                if(ss.getSize()<SS_MAX_SIZE){
                    ss.SSPush(pkt.first,updateValueToInt,0);
                }else{
                    double temp=(double)rand()/RAND_MAX;
                    uint32_t minVal=ss.getMinVal();
                    uint32_t newVal=minVal+updateValueToInt;

                    if(temp<(double)updateValueToInt/newVal){
                        ss.SSPush(pkt.first,newVal,minVal);
                    }
                    else{
                        SSKeyNode* minKeyNode=ss.getMinKeyNode();
                        ss.SSUpdate(minKeyNode,newVal);
                    }
                }
            }
        }
    }

    void getEstimatedFlowSpreads(unordered_map<uint32_t, uint32_t>& estimatedFlowSpreads) const{
        ss.getKeyVals(estimatedFlowSpreads);
    }

    ~SKETCH(){
        delete []regArray;
    }
};

bool cmpPairFunc(pair<uint32_t,uint32_t>p1, pair<uint32_t,uint32_t>p2)
{
    return p1.second > p2.second;
}

vector<vector<double>> calculateMetrics(unordered_map<uint32_t, uint32_t> &estFlowSpreads,unordered_map<uint32_t, uint32_t> &actualFlowSpreads,vector<uint32_t> topKs, string outputDirPath,string savedFileName,string info){
    vector<pair<uint32_t,uint32_t>> estFlowSpreadsVec;
    for(auto iter=estFlowSpreads.begin();iter!=estFlowSpreads.end();iter++){
        estFlowSpreadsVec.push_back(make_pair(iter->first,iter->second));
    }
    sort(estFlowSpreadsVec.begin(), estFlowSpreadsVec.end(), cmpPairFunc);

    vector<pair<uint32_t,uint32_t>> actualFlowSpreadsVec;
    for(auto iter=actualFlowSpreads.begin();iter!=actualFlowSpreads.end();iter++){
        actualFlowSpreadsVec.push_back(make_pair(iter->first,iter->second));
    }
    sort(actualFlowSpreadsVec.begin(), actualFlowSpreadsVec.end(), cmpPairFunc);

    string resultFilePath = outputDirPath+savedFileName + ".result";
    ofstream resultFile;
    resultFile.open(resultFilePath,ios::out);

    resultFile<<"flowId\t\tacutalSpread\t\testSpread\t\tsort by est"<<endl;

    for(uint32_t i=0;i<estFlowSpreadsVec.size();i++){
        uint32_t key=estFlowSpreadsVec[i].first;
        uint32_t estSpread=estFlowSpreadsVec[i].second;
        uint32_t actualSpread=1.0;

        auto iter=actualFlowSpreads.find(key);
        if(iter!=actualFlowSpreads.end()){
            actualSpread=actualFlowSpreads[key];
        }

        resultFile<<key<<"\t\t"<<actualSpread<<"\t\t"<<estSpread<<endl;
    }
    resultFile.close();

    string metricsFilePath = outputDirPath+savedFileName + ".metrics";
    ofstream metricsFile;
    metricsFile.open(metricsFilePath,ios::app);

    metricsFile<<"******************"<<endl;
    metricsFile<<"Super Spreader Metrics: "<<info<<endl;
    metricsFile<<"******************"<<endl;
    //trueFlowNum means the true number of super spreaders with that threshold
    //fakeFlowNum means the number of fake flows that returned by reversible calculation operation based on Chinese Remainder Theorem, in other words, these flows do not exist in the packet stream.
    metricsFile<<"Threshold\t\tPR\tRC\tF1\tARE\tAAE\tTP\tFP\tFN\ttrueFlowNum\t\tfakeFlowNum"<<endl;

    vector<double> THs,TPs,FPs,FNs,PRs,RCs,F1s,AREs,AAEs,realFlowNums,fakeFlowNums;
    double fakeFlowNum=0;

    for(uint32_t k:topKs){
        double threshold=0.5*(actualFlowSpreadsVec[k-1].second+actualFlowSpreadsVec[k].second);

        unordered_set<uint32_t> trueSuperSpreads;
        for(uint32_t i=0;i<actualFlowSpreadsVec.size();i++){
            if(actualFlowSpreadsVec[i].second>=threshold){
                trueSuperSpreads.insert(actualFlowSpreadsVec[i].first);
            }
        }

        double TP=0,FP=0,totalRE=0,totalAE=0;

        for(uint32_t i=0;i<estFlowSpreadsVec.size();i++){
            uint32_t key=estFlowSpreadsVec[i].first;
            double estSpread=estFlowSpreadsVec[i].second;

            if(estSpread>=threshold){
                auto iter=trueSuperSpreads.find(key);
                if(iter!=trueSuperSpreads.end()){
                    TP+=1;
                }else{
                    FP+=1;
                }

                double actualSpread=1.0;
                auto iter2=actualFlowSpreads.find(key);
                if(iter2!=actualFlowSpreads.end()){
                    actualSpread=actualFlowSpreads[key];
                }else{
                    fakeFlowNum+=1;
                }

                totalAE+=abs(estSpread-actualSpread);
                totalRE+=abs((estSpread-actualSpread)/actualSpread);
            }else{
                break;
            }
        }

        double FN=trueSuperSpreads.size()-TP;
        double PR=TP/(TP+FP);
        double RC=TP/(TP+FN);
        double F1=2*PR*RC/(PR+RC);
        double ARE=totalRE/(TP+FP);
        double AAE=totalAE/(TP+FP);

        char temp[500]{0};

        sprintf(temp,"%-10.3f\t%-10.3f\t%-10.3f\t%-10.3f\t%-10.3f\t%-10.3f\t%-6d\t%-6d\t%-6d\t%-6d\t\t%-6d",threshold,PR,RC,F1,ARE,AAE,int(TP),int(FP),int(FN),trueSuperSpreads.size(),int(fakeFlowNum));
        sprintf(temp,"threshold:%-10.3f\t PR:%-10.3f\t RC:%-10.3f\t F1:%-10.3f\t ARE:%-10.3f\t TP:%-6d\t FP:%-6d\t FN:%-6d\t",threshold,PR,RC,F1,ARE,int(TP),int(FP),int(FN));
//         sprintf(temp,"%-10.3f\t PR:%-10.3f\t RC:%-10.3f\t F1:%-10.3f\t ARE:%-10.3f\t AAE:%-10.3f\t %-6d\t%-6d\t%-6d\t%-6d\t\t%-6d",threshold,PR,RC,F1,ARE,AAE,int(TP),int(FP),int(FN),trueSuperSpreads.size(),int(fakeFlowNum));
        metricsFile<<temp<<endl;

        THs.push_back(threshold);
        PRs.push_back(PR);
        RCs.push_back(RC);
        F1s.push_back(F1);
        AREs.push_back(ARE);
        AAEs.push_back(AAE);
        TPs.push_back(TP);
        FPs.push_back(FP);
        FNs.push_back(FN);
        realFlowNums.push_back(trueSuperSpreads.size());
        fakeFlowNums.push_back(fakeFlowNum);
    }

    metricsFile.close();

    vector<vector<double>> metrics={THs,PRs,RCs,F1s,AREs,AAEs,TPs,FPs,FNs,realFlowNums,fakeFlowNums};
    return metrics;
}

//for single measurement point
unsigned int ReadInSingleTraceFile(const char* traceFilePath, vector<pair<uint32_t,uint32_t>>& pkts, unordered_map<uint32_t, uint32_t>& actualFlowSpreads,unordered_map<uint32_t, unordered_set<uint32_t>>& allFlowElements)
{

    FILE* fin = fopen(traceFilePath, "rb");

    uint64_t aPkt;
    unsigned int count = 0;

    while (fread(&aPkt, 8, 1, fin) == 1) {
        uint32_t src=aPkt&0xFFFFFFFF;
        uint32_t dst=aPkt>>32;

        pkts.push_back(make_pair(src,dst));
        allFlowElements[src].insert(dst);

        count++;
        if (count % 5000000 == 0) {
            printf("Successfully read in %s, %u packets\n", traceFilePath, count);

        }
    }
    printf("Successfully read in %s, %u packets\n", traceFilePath, count);
    fclose(fin);

    printf("Successfully get actual flow spreads\n");

    return count;
}

void getDataSets(const char* traceFileDir,uint32_t fileNum,vector<vector<pair<uint32_t,uint32_t>>>& pktsOfFiles,vector<unordered_map<uint32_t, uint32_t>>& actualFlowSpreadsOfFiles,vector<uint32_t>& actualItemNumOfFiles){
    unordered_map<uint32_t, uint32_t> actualFlowSpreads;
    unordered_map<uint32_t, unordered_set<uint32_t>> allFlowElements;
    for(uint32_t i=0;i<fileNum;i++) {
        char traceFilePath[256]{0};
        sprintf(traceFilePath, "%s%02d.dat", traceFileDir, i);

        //prepare dataset
        cout << "prepare dataset: " << traceFilePath << endl;
        vector<pair<uint32_t, uint32_t>> pkts;
        unsigned int actualItemNum = ReadInSingleTraceFile(traceFilePath, pkts, actualFlowSpreads,allFlowElements);

        cout << "number of packets: " << actualItemNum << endl;
        // cout << "number of flows: " << actualFlowSpreads.size() << endl;
        cout << "*********************" << endl;

        pktsOfFiles.push_back(pkts);
        // actualFlowSpreadsOfFiles.push_back(actualFlowSpreads);
        actualItemNumOfFiles.push_back(actualItemNum);
    }
    for(auto iter=allFlowElements.begin();iter!=allFlowElements.end();iter++){
        uint32_t key=iter->first;
        uint32_t actualSpread=iter->second.size();
        if(actualFlowSpreads[key]==0){
            actualFlowSpreads[key]=actualSpread;
        }
        else{
            actualFlowSpreads[key]+=actualSpread;
        }
    }
    actualFlowSpreadsOfFiles.push_back(actualFlowSpreads);
    cout<<"finished"<<endl;
}

void processMultipleTracesAndGetAvgMetrics(uint32_t SS_MAX_SIZE,uint32_t REG_NUM,const char* traceFileDir,uint32_t fileNum,vector<vector<pair<uint32_t,uint32_t>>>& pktsOfFiles,vector<unordered_map<uint32_t, uint32_t>>& actualFlowSpreadsOfFiles,vector<uint32_t>& actualItemNumOfFiles,const vector<uint32_t>& topKs, string outputDirPath,string traceInfo){
    //in the stream summary, suppose each pointer is 32 bits
    //each entry of the hash table only stores a pointer (SSKeyNode*): 32 bits
    //each key node stores a key, an error counter, 2 pointers (one to the parent value node, one to the next key node): 32+32+32*2=128 bits
    //each value node stores a value, 3 pointers (to the first key node, the pre value node, and the next value node): 32+32*3=128 bits
    //suppose each flow needs one key node, one value node, and one hash table entry (simply suppose the fill factor of the hash table is 1), then each flow needs 32+128+128= 288 bits

    uint32_t totalMem=SS_MAX_SIZE*288+REG_NUM*5;
    double SSMemRatio=(double)SS_MAX_SIZE*288/totalMem;
    cout << "totalMem:" << totalMem/(8*1024.0) << "KB" << endl;
    cout << "*********************" << endl;

    char temp[500]{0};
    sprintf(temp,"%s data %d files src FreeRS_INC_DE SSSz=%d RNum=%d SSR=%.2f mem=%.2fKB",traceInfo.c_str(),fileNum,SS_MAX_SIZE,REG_NUM,SSMemRatio,totalMem/(8.0*1024));
    string savedFileName(temp);
    
    double avgThroughput=0;
    vector<double> THs(topKs.size()),TPs(topKs.size()),FPs(topKs.size()),FNs(topKs.size()),PRs(topKs.size());
    vector<double> RCs(topKs.size()),F1s(topKs.size()),AREs(topKs.size()),AAEs(topKs.size()),realFlowNums(topKs.size()),fakeFlowNums(topKs.size());
    vector<vector<double>> avgMetrics={THs,PRs,RCs,F1s,AREs,AAEs,TPs,FPs,FNs,realFlowNums,fakeFlowNums};

    char traceFilePath[256]{0};
    sprintf(traceFilePath,"%s%02d--%02d.dat",traceFileDir,0,fileNum-1);
    
    SKETCH sketch(SS_MAX_SIZE,REG_NUM);
    
    unordered_map<uint32_t, uint32_t> actualFlowSpreads=actualFlowSpreadsOfFiles[0];
    unsigned int actualItemNum = 0;
    cout<<"fi 2"<<endl;
    string metricsFilePath = outputDirPath+savedFileName + "-avg.metrics";
    ofstream metricsFile;
    clock_t time1 = clock();
    
    cout<<"finish 4"<<endl;
    
    for(uint32_t i=0;i<fileNum;i++){

        vector<pair<uint32_t,uint32_t>> pkts=pktsOfFiles[i];
        unsigned int actualItemNumPart = actualItemNumOfFiles[i];
        cout<<"pktsSize"<<pkts.size()<<endl;
        cout << "prepare algorithm "<< endl;
        for (unsigned int i = 0; i < pkts.size(); i++) {
            sketch.insert(pkts[i]);
        }
        actualItemNum+=actualItemNumPart;
    }
    clock_t time2 = clock();
    double numOfSeconds = (double)(time2 - time1) / CLOCKS_PER_SEC;//the seconds using to insert items
    double throughput = (actualItemNum / 1000000.0) / numOfSeconds;
    cout << "throughput: " << throughput << " Mpps" << ", each insert operation uses " << 1000.0 / throughput << " ns" << endl;
    cout << "*********************" << endl;

    unordered_map<uint32_t, uint32_t> estFlowSpreads;
    sketch.getEstimatedFlowSpreads(estFlowSpreads);

    vector<vector<double>> metrics=calculateMetrics(estFlowSpreads,actualFlowSpreads,topKs,outputDirPath,savedFileName+"-avg",traceFilePath);
    metricsFile.open(metricsFilePath,ios::app);
    metricsFile<<"throughput: " << throughput << " Mpps" << ", each insert operation uses " << 1000.0 / throughput << " ns" << endl;
    metricsFile.close();

    avgThroughput+=throughput;
    for(uint32_t metricIdx=0;metricIdx<metrics.size();metricIdx++){
        for(uint32_t kIdx=0;kIdx<topKs.size();kIdx++){
            avgMetrics[metricIdx][kIdx]+=metrics[metricIdx][kIdx];
        }
    }
    
    metricsFile.open(metricsFilePath,ios::app);

    metricsFile<<"******************"<<endl;
    metricsFile<<"Average Super Spreader Metrics"<<endl;
    metricsFile<<"******************"<<endl;
    //trueFlowNum means the true number of super spreaders with that threshold
    //fakeFlowNum means the number of fake flows that returned by reversible calculation operation based on Chinese Remainder Theorem, in other words, these flows do not exist in the packet stream.
    metricsFile<<"Threshold\t\tPR\tRC\tF1\tARE\tAAE\tTP\tFP\tFN\ttrueFlowNum\t\tfakeFlowNum"<<endl;

    //THs,PRs,RCs,F1s,AREs,AAEs,TPs,FPs,FNs,realFlowNums,fakeFlowNums};
    for(uint32_t kIdx=0;kIdx<topKs.size();kIdx++){
        for(uint32_t metricIdx=0;metricIdx<avgMetrics.size();metricIdx++){

            avgMetrics[metricIdx][kIdx]=avgMetrics[metricIdx][kIdx]/fileNum;
            char temp[500]{0};
            sprintf(temp,"%-10.3f\t",avgMetrics[metricIdx][kIdx]);
            metricsFile<<temp;
        }
        metricsFile<<endl;

        if(kIdx==topKs.size()-1){
            avgF1RatioFile.open(avgF1RatioFileName,ios::app);
            char temp[50]{0};
            sprintf(temp,"%.3f",avgMetrics[3][kIdx]);
            avgF1RatioFile<<temp<<",";
            avgF1RatioFile.close();
        }
    }
    avgThroughput=avgThroughput/fileNum;
    metricsFile<<"Average throughput: " << avgThroughput << " Mpps" << endl;
    metricsFile.close();

    cout << "Average throughput: " << avgThroughput << " Mpps" << endl;
    cout <<metricsFilePath<<endl;
}



int main() {

    vector<uint32_t> topKs={100,150,200,250,300};
    string outputDirPath="./resultsFreeRS/2019/";//test
//     string outputDirPath="./resultsFreeRS/2019/";//test
//     string outputDirPath="./resultsFreeRS/Cos/";//test
//     string outputDirPath="./resultsFreeRS/20204/";//test  
//     string outputDirPath="./resultsFreeRS/Mul/";//test 
    
    //     const char* traceFileDir=R"(../data/2016ipv4/)";
    const char* traceFileDir=R"(../data/2019ipv4/)";
//     const char* traceFileDir=R"(../data/Mul/)";
//     const char* traceFileDir=R"(/home/tmp/Desktop/ZBY/throughput/perSource_LittleEndian/data/Mul2020.4/)";
//     const char* traceFileDir=R"(/home/tmp/Desktop/ZBY/throughput/perSource_LittleEndian/data/Cos/)";

    uint32_t fileNum=10;
    vector<vector<pair<uint32_t,uint32_t>>> pktsOfFiles;
    vector<unordered_map<uint32_t, uint32_t>> actualFlowSpreadsOfFiles;
    vector<uint32_t> actualItemNumOfFiles;

//    getDataSets(traceFileDir2016,fileNum,pktsOfFiles2016,actualFlowSpreadsOfFiles2016,actualItemNumOfFiles2016);
//    getDataSets(traceFileDir2019,fileNum,pktsOfFiles2019,actualFlowSpreadsOfFiles2019,actualItemNumOfFiles2019);
//    getDataSets(traceFileDirCos,fileNum,pktsOfFilesCos,actualFlowSpreadsOfFilesCos,actualItemNumOfFilesCos);
    getDataSets(traceFileDir,fileNum,pktsOfFiles,actualFlowSpreadsOfFiles,actualItemNumOfFiles);

    for(int mem=100;mem<=300;mem+=50){

        avgF1RatioFile.open(avgF1RatioFileName,ios::app);
        avgF1RatioFile<<"mem="<<mem<<" KB"<<endl;
        avgF1RatioFile.close();

        for(int i=80;i<=80;i+=2){
            double SSMemRatio=i*0.01;
            uint32_t totalMemBits=mem*8*1024;
            double SSMem=totalMemBits*SSMemRatio;
            uint32_t SSSize=uint32_t(SSMem/288.0);
            double regMem=totalMemBits-SSSize*288;
            uint32_t regNum=uint32_t(regMem/5.0);


//            processMultipleTracesAndGetAvgMetrics(SSSize,regNum,traceFileDir2016,fileNum,pktsOfFiles2016,actualFlowSpreadsOfFiles2016,actualItemNumOfFiles2016,topKs,outputDirPath,"2016");
//            processMultipleTracesAndGetAvgMetrics(SSSize,regNum,traceFileDir2019,fileNum,pktsOfFiles2019,actualFlowSpreadsOfFiles2019,actualItemNumOfFiles2019,topKs,outputDirPath,"2019");
            processMultipleTracesAndGetAvgMetrics(SSSize,regNum,traceFileDir,fileNum,pktsOfFiles,actualFlowSpreadsOfFiles,actualItemNumOfFiles,topKs,outputDirPath,"2019");
//            processMultipleTracesAndGetAvgMetrics(SSSize,regNum,traceFileDirMul,fileNum,pktsOfFilesMul,actualFlowSpreadsOfFilesMul,actualItemNumOfFilesMul,topKs,outputDirPath,"Mul_2020_4");

        }
        avgF1RatioFile.open(avgF1RatioFileName,ios::app);
        avgF1RatioFile<<endl;
        avgF1RatioFile.close();
    }

    return 0;
}
