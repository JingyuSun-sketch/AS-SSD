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
#include <numeric>
#include <fstream>
#include <string>
#include "MurmurHash3.h"

using namespace std;
const uint32_t HASH_SEEDS10[]={2725628290,1444492164,1001369998,2291192083,3632694421,3217476757,962572696,3523028760,1564808096,1686744101};
const uint32_t HASH_SEEDS[]={674876105,860938402,3044024791,706198233,1838710314,2602646517};

class BITMAP{
private:
    uint8_t* bitmap;
    uint32_t size;
    uint32_t oneBitsNum=0;
public:
    BITMAP(uint32_t size){
        uint32_t byteNum=ceil(size/8.0);
        bitmap=new uint8_t[byteNum]();
        this->size=size;
    }
    BITMAP(const BITMAP& c){
        uint32_t byteNum=ceil(c.size/8.0);
        bitmap=new uint8_t[byteNum]();
        size=c.size;
        for(uint32_t byteIdx=0;byteIdx<byteNum;byteIdx++){
            bitmap[byteIdx]=c.bitmap[byteIdx];
        }
        oneBitsNum=c.oneBitsNum;
    }
    BITMAP& operator=(const BITMAP& other){
        if(this==&other){
            return *this;
        }
        delete[] bitmap;

        uint32_t byteNum=ceil(other.size/8.0);
        bitmap=new uint8_t[byteNum]();
        size=other.size;
        for(uint32_t byteIdx=0;byteIdx<byteNum;byteIdx++){
            bitmap[byteIdx]=other.bitmap[byteIdx];
        }
        oneBitsNum=other.oneBitsNum;
        return *this;
    }

    bool test(uint32_t bitIdx) const{
        uint32_t byteIdx=bitIdx/8;
        uint32_t bitIdxInByte=bitIdx%8;
        uint8_t mask=1<<bitIdxInByte;

        if((bitmap[byteIdx]&mask)!=0){
            return true;
        }else{
            return false;
        }
    }

    void set(uint32_t bitIdx){
        uint32_t byteIdx=bitIdx/8;
        uint32_t bitIdxInByte=bitIdx%8;
        uint8_t mask=1<<bitIdxInByte;

        if((bitmap[byteIdx]&mask)==0){
            bitmap[byteIdx]=bitmap[byteIdx]|mask;
            oneBitsNum++;
        }
    }
    void reset(uint32_t bitIdx){
        uint32_t byteIdx=bitIdx/8;
        uint32_t bitIdxInByte=bitIdx%8;
        uint8_t mask=1<<bitIdxInByte;

        if((bitmap[byteIdx]&mask)!=0){
            bitmap[byteIdx]=bitmap[byteIdx]&(~mask);
            oneBitsNum--;
        }
    }

    void intersectOp(const BITMAP& other){
        for(uint32_t bitIdx=0;bitIdx<size;bitIdx++){
            if(!other.test(bitIdx) && test(bitIdx)){
                reset(bitIdx);
            }
        }
    }
    void unionOp(const BITMAP& other){
        for(uint32_t bitIdx=0;bitIdx<size;bitIdx++){
            if(other.test(bitIdx) && !test(bitIdx)){
                set(bitIdx);
            }
        }
    }

    uint32_t getSize() const{
        return size;
    }
    double getOneBitsRatio() const{
        return (double)oneBitsNum/size;
    }
    double getZeroBitsRatio() const{
        return (double)(size-oneBitsNum)/size;
    }

    ~BITMAP(){
        delete[] bitmap;
    }
};

template<uint32_t BASE_BIT_NUM,uint8_t PIMAX>
class BUCKET{
public:
    BITMAP EC;
    uint8_t ET=0;
    uint32_t PS=0;
    uint8_t PI=0;
    vector<unordered_set<uint32_t>> QS;
    BUCKET():EC(BASE_BIT_NUM),QS(BASE_BIT_NUM){};
    void update(const pair<uint32_t,uint32_t>& pkt, double delta,uint32_t HASH_SEED){
        uint32_t hashValue1=0;
        MurmurHash3_x86_32(&pkt,8,HASH_SEED, &hashValue1);
        uint32_t bitNum=BASE_BIT_NUM<<ET;
        uint32_t quotient=hashValue1/bitNum;
        uint32_t bitIdx=hashValue1%bitNum;
        EC.set(bitIdx);

        QS[bitIdx].insert(quotient);

        if(PI==0){
            PS=pkt.first;
            PI++;
        }else{
            if(PS==pkt.first){
                if(PI<PIMAX){
                    PI++;
                }
            }else{
                double randVal=(double)rand()/RAND_MAX;
                if(randVal<delta){
                    PI--;
                }
            }

        }
    }
    void extend(){
        uint32_t curBitNum=BASE_BIT_NUM<<ET;
        uint32_t newBitNum=BASE_BIT_NUM<<(ET+1);

        BITMAP newEC(newBitNum);
        vector<unordered_set<uint32_t>> newQS(newBitNum);

        for(uint32_t bitIdx=0;bitIdx<curBitNum;bitIdx++){
            for(auto iter=QS[bitIdx].begin();iter!=QS[bitIdx].end();iter++){
                uint32_t quotient=*iter;
                if(quotient%2==0){
                    newEC.set(bitIdx);
                    newQS[bitIdx].insert(quotient/2);
                }else{
                    newEC.set(bitIdx+curBitNum);
                    newQS[bitIdx+curBitNum].insert((quotient-1)/2);
                }
            }
        }

        EC=newEC;
        ET++;
        QS=newQS;
    }

    double getOneBitsRatio() const{
        return EC.getOneBitsRatio();
    }

    uint32_t query() const{
        uint32_t curBitNum=BASE_BIT_NUM<<ET;
        double emptyRatio=EC.getZeroBitsRatio();

        if(emptyRatio==0){
            emptyRatio=1.0/curBitNum;
        }
        double val=-1.0*curBitNum*log(emptyRatio);
        return round(val);
    }
};

template<uint32_t ROW_NUM,uint32_t COL_NUM, uint32_t BASE_BIT_NUM,uint8_t PIMAX>
class SKETCH{
private:
    vector<vector<BUCKET<BASE_BIT_NUM,PIMAX>>> buckets;
    double epsilon;
    double delta;

public:
    SKETCH(double epsilon,double delta):buckets(ROW_NUM,vector<BUCKET<BASE_BIT_NUM,PIMAX>>(COL_NUM)){
        this->epsilon=epsilon;
        this->delta=delta;
    }

    void insert(const pair<uint32_t,uint32_t>& pkt,uint32_t HASH_SEED) {
        uint32_t colIdxes[ROW_NUM]{0};
        uint32_t fullBucketNum=0;
        for (uint32_t rowIdx = 0; rowIdx < ROW_NUM; rowIdx++) {
            uint32_t hashValue = 0;
            MurmurHash3_x86_32(&pkt.first, 4, HASH_SEEDS[rowIdx], &hashValue);
            uint32_t colIdx = hashValue % COL_NUM;
            colIdxes[rowIdx]=colIdx;

            if(buckets[rowIdx][colIdx].getOneBitsRatio()>=epsilon){
                fullBucketNum++;
            }
        }

        for (uint32_t rowIdx = 0; rowIdx < ROW_NUM; rowIdx++) {
            uint32_t colIdx=colIdxes[rowIdx];
            BUCKET<BASE_BIT_NUM,PIMAX>& bucket=buckets[rowIdx][colIdx];

            if(fullBucketNum>=ROW_NUM){
                bucket.extend();
            }

            bucket.update(pkt,delta,HASH_SEED);
        }
    }

//    //the estimation method given by the paper
//    //it leads to serious over-estimation.
//    //may be the paper write wrong query methods
//    uint32_t query(uint32_t key) const{
//        uint32_t minVal=-1;
//        for (uint32_t rowIdx = 0; rowIdx < ROW_NUM; rowIdx++) {
//            uint32_t hashValue = 0;
//            MurmurHash3_x86_32(&key, 4, HASH_SEEDS[rowIdx], &hashValue);
//            uint32_t colIdx = hashValue % COL_NUM;
//            uint32_t val=buckets[rowIdx][colIdx].query();
//            if(val<minVal){
//                minVal=val;
//            }
//        }
//        return minVal;
//    }

    // use bitwise-and operation like extendedSketch
    uint32_t query(uint32_t key) const{
        uint32_t hashValue = 0;
        uint32_t colIdx=0;
        uint32_t maxExtendTimes=0;

        uint32_t colIdxes[ROW_NUM]{0};
        for (uint32_t rowIdx = 0; rowIdx < ROW_NUM; rowIdx++) {
            MurmurHash3_x86_32(&key, 4, HASH_SEEDS[rowIdx], &hashValue);
            colIdx = hashValue % COL_NUM;
            colIdxes[rowIdx]=colIdx;
            BUCKET<BASE_BIT_NUM, PIMAX> bucket = buckets[rowIdx][colIdx];

            if(bucket.ET>maxExtendTimes){
                maxExtendTimes=bucket.ET;
            }
        }

        BUCKET<BASE_BIT_NUM,PIMAX> tempBucket=buckets[0][colIdxes[0]];
        while(tempBucket.ET<maxExtendTimes){
            tempBucket.extend();
        }

        for (uint32_t rowIdx = 1; rowIdx < ROW_NUM; rowIdx++) {


            BUCKET<BASE_BIT_NUM, PIMAX> bucket = buckets[rowIdx][colIdxes[rowIdx]];

            while(bucket.ET < maxExtendTimes)
            {
                bucket.extend();
            }

            tempBucket.EC.intersectOp(bucket.EC);
        }
        uint32_t val=tempBucket.query();
        return val;
    }

    void getEstimatedFlowSpreads(unordered_map<uint32_t, uint32_t>& estimatedFlowSpreads) const{
        for(uint32_t rowIdx=0;rowIdx<ROW_NUM;rowIdx++){
            for(uint32_t colIdx=0;colIdx<COL_NUM;colIdx++){
                uint32_t key=buckets[rowIdx][colIdx].PS;
                if(estimatedFlowSpreads.find(key)==estimatedFlowSpreads.end()){
                    uint32_t val= query(key);
                    estimatedFlowSpreads[key]=val;
                }
            }
        }
    }

    static void getMergedEstimatedFlowSpreads(const vector<SKETCH>& sketches,unordered_map<uint32_t, uint32_t>& estimatedFlowSpreads){
        SKETCH<ROW_NUM,COL_NUM,BASE_BIT_NUM,PIMAX> mergedSketch=sketches[0];

        for (uint32_t rowIdx = 0; rowIdx < ROW_NUM; rowIdx++) {
            for (uint32_t colIdx = 0; colIdx < COL_NUM; colIdx++) {
                    BUCKET<BASE_BIT_NUM,PIMAX>& mergedBucket=mergedSketch.buckets[rowIdx][colIdx];
                for(uint32_t sketchIdx=1;sketchIdx<sketches.size();sketchIdx++){
                    BUCKET<BASE_BIT_NUM,PIMAX> bucket=sketches[sketchIdx].buckets[rowIdx][colIdx];

                    while(mergedBucket.ET<bucket.ET){
                        mergedBucket.extend();
                    }
                    while(bucket.ET<mergedBucket.ET){
                        bucket.extend();
                    }
                    mergedBucket.EC.unionOp(bucket.EC);
                    uint32_t curBitNum=BASE_BIT_NUM<<bucket.ET;
                    for(uint32_t bitIdx=0;bitIdx<curBitNum;bitIdx++){
                        for(auto iter=bucket.QS[bitIdx].begin();iter!=bucket.QS[bitIdx].end();iter++){
                            mergedBucket.QS[bitIdx].insert(*iter);
                        }
                    }

                    if(bucket.PI>mergedBucket.PI){
                        mergedBucket.PI=bucket.PI;
                        mergedBucket.PS=bucket.PS;
                    }
                }
            }
        }

        mergedSketch.getEstimatedFlowSpreads(estimatedFlowSpreads);
    }

};

bool cmpPairFunc(pair<uint32_t,uint32_t>p1, pair<uint32_t,uint32_t>p2)
{
    return p1.second > p2.second;
}

vector<vector<double>> calculateMetrics(unordered_map<uint32_t, uint32_t> &estFlowSpreads,unordered_map<uint32_t, uint32_t> &actualFlowSpreads,vector<uint32_t> topKs, string outputDirPath,string savedFileName,string info,unordered_map<double, vector<double>>& avgF1){
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
        
        avgF1[k].push_back(F1);
        
        double ARE=totalRE/(TP+FP);
        double AAE=totalAE/(TP+FP);

        char temp[500]{0};
        sprintf(temp,"%-10.3f\t%-10.3f\t%-10.3f\t%-10.3f\t%-10.3f\t%-10.3f\t%-6d\t%-6d\t%-6d\t%-6d\t\t%-6d",threshold,PR,RC,F1,ARE,AAE,int(TP),int(FP),int(FN),trueSuperSpreads.size(),int(fakeFlowNum));
        sprintf(temp,"threshold:%-10.3f\t PR:%-10.3f\t RC:%-10.3f\t F1:%-10.3f\t ARE:%-10.3f\t TP:%-6d\t FP:%-6d\t FN:%-6d\t",threshold,PR,RC,F1,ARE,int(TP),int(FP),int(FN));
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

void getDataSets(const char* traceFileDir,uint32_t fileNum,uint32_t step,vector<vector<pair<uint32_t,uint32_t>>>& pktsOfFiles,vector<unordered_map<uint32_t, uint32_t>>& actualFlowSpreadsOfFiles,vector<uint32_t>& actualItemNumOfFiles){
    uint32_t fileNumMax=fileNum+1;
    for(uint32_t i=step;i<fileNumMax;i=i+step) {
        cout<<"i="<<i<<endl;
        unordered_map<uint32_t, uint32_t> actualFlowSpreads;
        unordered_map<uint32_t, unordered_set<uint32_t>> allFlowElements;
        for(uint32_t j=(i-step);j<i;j++){
            char traceFilePath[256]{0};
            sprintf(traceFilePath, "%s%02d.dat", traceFileDir, j);

            //prepare dataset
            cout << "prepare dataset: " << traceFilePath << endl;
            vector<pair<uint32_t, uint32_t>> pkts;
            
            unsigned int actualItemNum = ReadInSingleTraceFile(traceFilePath, pkts, actualFlowSpreads,allFlowElements);

            cout << "number of packets: " << actualItemNum << endl;
            cout << "*********************" << endl;

            pktsOfFiles.push_back(pkts);
            // cout<<" pktsOfFiles:"<<pktsOfFiles.size()<<endl;
            
            // cout<<" actualFlowSpreadsOfFiles:"<<actualFlowSpreadsOfFiles.size()<<endl;
            actualItemNumOfFiles.push_back(actualItemNum);
        }
        for(auto iter=allFlowElements.begin();iter!=allFlowElements.end();iter++){
            uint32_t key=iter->first;
            uint32_t acutalSpread=iter->second.size();
            if(actualFlowSpreads[key]==0){
                actualFlowSpreads[key]=acutalSpread;
            }
            else{
                actualFlowSpreads[key]+=acutalSpread;
            }
        }
        cout << "number of flows: " << actualFlowSpreads.size() << endl;
        cout << "*********************" << endl;
        actualFlowSpreadsOfFiles.push_back(actualFlowSpreads);
        cout<<"one turn finished,and i="<<i<<endl;
    }
    cout<<"finished"<<endl;
}

template<uint32_t ROW_NUM,uint32_t COL_NUM, uint32_t BASE_BIT_NUM,uint8_t PIMAX>
void processMultipleTracesAndGetAvgMetrics(const char* traceFileDir,uint32_t fileNum,uint32_t step,vector<vector<pair<uint32_t,uint32_t>>>& pktsOfFiles,vector<unordered_map<uint32_t, uint32_t>>& actualFlowSpreadsOfFiles,vector<uint32_t>& actualItemNumOfFiles,const vector<uint32_t>& topKs, string outputDirPath,string traceInfo,double epsilon,double delta){
    //PS: 32bits, PI: log2(PIMAX+1) bits, EC bitmap: baseBitsNum bits
    //recording the number of used bits in bitmap and the number of extended times use 32 bits and 8 bits separately
    //the pointer pointing to the extensible bitmap uses 32bits
    //each pointer pointing to the quotient set for each bit in the bitmap uses 32 bits
    uint32_t totalMem=ROW_NUM*COL_NUM*(32+ceil(log2(PIMAX+1))+32+BASE_BIT_NUM*(1+32)+32+8);
    cout << "totalMem:" << totalMem/(8*1024.0) << "KB" << endl;
    cout << "*********************" << endl;

    char temp[500]{0};
    sprintf(temp,"%s data %d files src ESP r=%d w=%d b=%d ep=%.2f dt=%.2f PIMAX=%d mem=%.2fKB",traceInfo.c_str(),fileNum,ROW_NUM,COL_NUM,BASE_BIT_NUM,epsilon,delta,PIMAX,totalMem/(8.0*1024));
    string savedFileName(temp);

    double avgThroughput=0;
    vector<double> THs(topKs.size()),TPs(topKs.size()),FPs(topKs.size()),FNs(topKs.size()),PRs(topKs.size());
    vector<double> RCs(topKs.size()),F1s(topKs.size()),AREs(topKs.size()),AAEs(topKs.size()),realFlowNums(topKs.size()),fakeFlowNums(topKs.size());
    vector<vector<double>> avgMetrics={THs,PRs,RCs,F1s,AREs,AAEs,TPs,FPs,FNs,realFlowNums,fakeFlowNums};
    string metricsFilePath = outputDirPath+savedFileName + "-avg.metrics";
    ofstream metricsFile;
    unordered_map<double, vector<double>> avgF1;
    for(int eachTime=0;eachTime<1;eachTime++){
        uint32_t HASH_SEED = HASH_SEEDS10[eachTime];
        cout<<HASH_SEED<<endl;
        for(uint32_t i=step;i<=fileNum;i=i+step){
            char traceFilePath[256]{0};
            sprintf(traceFilePath,"%s%02d.dat--%02d.dat",traceFileDir,i-step,i-1);
            unordered_map<uint32_t, uint32_t> actualFlowSpreads=actualFlowSpreadsOfFiles[(i/step)-1];
            unsigned int actualItemNum=0;

        //prepare dataset
            cout << "number of packets: " << actualItemNum << endl;
            cout << "number of flows: " << actualFlowSpreads.size() << endl;
            cout << "*********************" << endl;


            cout << "prepare algorithm "<< endl;
            SKETCH<ROW_NUM,COL_NUM,BASE_BIT_NUM,PIMAX> sketch(epsilon,delta);

            clock_t time1 = clock();
            for(uint32_t j=(i-step);j<i;j++){
                vector<pair<uint32_t,uint32_t>> pkts=pktsOfFiles[j];
                unsigned int actualItemNumPart = actualItemNumOfFiles[j];
                cout<<"pktsSize"<<pkts.size()<<endl;
                cout << "prepare algorithm "<< endl;       
                for (unsigned int i = 0; i < pkts.size(); i++) {
                    sketch.insert(pkts[i],HASH_SEED);
                }
                actualItemNum+=actualItemNumPart;
            }       
            clock_t time2 = clock();
            cout<<"finished"<<endl;

            double numOfSeconds = (double)(time2 - time1) / CLOCKS_PER_SEC;//the seconds using to insert items
            double throughput = (actualItemNum / 1000000.0) / numOfSeconds;
            cout << "throughput: " << throughput << " Mpps" << ", each insert operation uses " << 1000.0 / throughput << " ns" << endl;
            cout << "*********************" << endl;
            unordered_map<uint32_t, uint32_t> estFlowSpreads;
            sketch.getEstimatedFlowSpreads(estFlowSpreads);

            vector<vector<double>> metrics=calculateMetrics(estFlowSpreads,actualFlowSpreads,topKs,outputDirPath,savedFileName+"-avg",traceFilePath,avgF1);

            metricsFile.open(metricsFilePath,ios::app);
            metricsFile<<"throughput: " << throughput << " Mpps" << ", each insert operation uses " << 1000.0 / throughput << " ns" << endl;
            metricsFile.close();

            avgThroughput+=throughput;
            for(uint32_t metricIdx=0;metricIdx<metrics.size();metricIdx++){
                for(uint32_t kIdx=0;kIdx<topKs.size();kIdx++){
                    avgMetrics[metricIdx][kIdx]+=metrics[metricIdx][kIdx];
                }
            }
        }
    }
    metricsFile.open(metricsFilePath,ios::app);
    for(auto iter=avgF1.begin();iter!=avgF1.end();iter++){
        double key=iter->first;
        vector<double> keyF1=iter->second;
        metricsFile<<"key:"<<key<<" avg_F1"<<accumulate(begin(keyF1),end(keyF1),0.0)/keyF1.size()<<endl;
    }
    metricsFile.close();
}


int main() {

    vector<uint32_t> topKs={100,200,300,400,500};

    string outputDirPath="./results/Five/";

    double epsilon=0.75;
    double delta=0.1;
    const uint32_t PIMAX=15;
    const uint32_t BASE_BIT_NUM=32;
    const uint32_t ROW_NUM=4;

//    //100KB
   const uint32_t COL_NUM_100KB=175;
//    //200KB
   const uint32_t COL_NUM_200KB=351;
//    //300KB
   const uint32_t COL_NUM_300KB=527;
//    //150KB
   const uint32_t COL_NUM_150KB=264;
//    //250KB
   const uint32_t COL_NUM_250KB=440;


    uint32_t fileNum=5;
    uint32_t step=5;
    const char* traceFileDir2019=R"(../data/2019ipv4/)";
    vector<vector<pair<uint32_t,uint32_t>>> pktsOfFiles;
    vector<unordered_map<uint32_t, uint32_t>> actualFlowSpreadsOfFiles;
    vector<uint32_t> actualItemNumOfFiles;
    getDataSets(traceFileDir2019,fileNum,step,pktsOfFiles,actualFlowSpreadsOfFiles,actualItemNumOfFiles);
    cout<<"get data set finished"<<endl;

   //single measurement point

   processMultipleTracesAndGetAvgMetrics<ROW_NUM,COL_NUM_100KB,BASE_BIT_NUM,PIMAX>(traceFileDir2019,fileNum,step,pktsOfFiles,actualFlowSpreadsOfFiles,actualItemNumOfFiles,topKs,outputDirPath,"2019",epsilon,delta);
   processMultipleTracesAndGetAvgMetrics<ROW_NUM,COL_NUM_200KB,BASE_BIT_NUM,PIMAX>(traceFileDir2019,fileNum,step,pktsOfFiles,actualFlowSpreadsOfFiles,actualItemNumOfFiles,topKs,outputDirPath,"2019",epsilon,delta);
   processMultipleTracesAndGetAvgMetrics<ROW_NUM,COL_NUM_300KB,BASE_BIT_NUM,PIMAX>(traceFileDir2019,fileNum,step,pktsOfFiles,actualFlowSpreadsOfFiles,actualItemNumOfFiles,topKs,outputDirPath,"2019",epsilon,delta);
   processMultipleTracesAndGetAvgMetrics<ROW_NUM,COL_NUM_150KB,BASE_BIT_NUM,PIMAX>(traceFileDir2019,fileNum,step,pktsOfFiles,actualFlowSpreadsOfFiles,actualItemNumOfFiles,topKs,outputDirPath,"2019",epsilon,delta);
   processMultipleTracesAndGetAvgMetrics<ROW_NUM,COL_NUM_250KB,BASE_BIT_NUM,PIMAX>(traceFileDir2019,fileNum,step,pktsOfFiles,actualFlowSpreadsOfFiles,actualItemNumOfFiles,topKs,outputDirPath,"2019",epsilon,delta);


    return 0;
}

