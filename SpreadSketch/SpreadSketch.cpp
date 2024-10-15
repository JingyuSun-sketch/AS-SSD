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
#include <bitset>

// #define HASH_SEED 51318471
//#define HASH_SEED 74151457
using namespace std;


const uint32_t HASH_SEEDS[]={674876105,860938402,3044024791,706198233,1838710314,2602646517};
const uint32_t HASH_SEEDS10[]={2725628290,1444492164,1001369998,2291192083,3632694421,3217476757,962572696,3523028760,1564808096,1686744101};
template<uint32_t BIT_NUM1,uint32_t BIT_NUM2,uint32_t SET_MAX,uint32_t SET_LAST_MAX,uint32_t BITMAP_NUM>
class MultiResolutionBitmap{
private:
    vector<bitset<BIT_NUM1>> baseBitmaps;
    bitset<BIT_NUM2> lastBitmap;
public:
    MultiResolutionBitmap():baseBitmaps(BITMAP_NUM-1){};

//    void update(uint32_t element){
//        uint32_t hashValue1=0;
//        MurmurHash3_x86_32(&element,4,HASH_SEED, &hashValue1);
//
//        uint32_t bitmapIdx=__builtin_clz(hashValue1);
//
//        uint32_t hashValue2=0;
//        MurmurHash3_x86_32(&element, 4,HASH_SEED+15417891, &hashValue2);
//        if(bitmapIdx<BITMAP_NUM-1){
//            uint32_t bitIdx=hashValue2%BIT_NUM1;
//            baseBitmaps[bitmapIdx].set(bitIdx);
//        }else{
//            uint32_t bitIdx=hashValue2%BIT_NUM2;
//            lastBitmap.set(bitIdx);
//        }
//    }
    void update(const pair<uint32_t,uint32_t>& pkt,uint32_t mszNum,uint32_t HASH_SEED){
        uint32_t hashValue2=0;
        MurmurHash3_x86_32(&pkt, 8,HASH_SEED+15417891, &hashValue2);
        if(mszNum<BITMAP_NUM-1){
            uint32_t bitIdx=hashValue2%BIT_NUM1;
            baseBitmaps[mszNum].set(bitIdx);
        }else{
            uint32_t bitIdx=hashValue2%BIT_NUM2;
            lastBitmap.set(bitIdx);
        }
    }
    uint32_t query() const{
        int32_t tempIdx=BITMAP_NUM-2;
        while(tempIdx>=0 && baseBitmaps[tempIdx].count()<=SET_MAX){
            tempIdx--;
        }
        uint32_t baseIdx=tempIdx+1;

        if(baseIdx==BITMAP_NUM-1 && lastBitmap.count()>SET_LAST_MAX){
            if(lastBitmap.count()>=BIT_NUM2){
                //the last bitmap is full, and cannot get estimated value by linear counting
                //return the estimation when there is 1 zero-bit
                double baseEst=BIT_NUM2*log(BIT_NUM2);
                double factor=pow(2,baseIdx);
                return round(baseEst*factor);
            }else{
                //not accurate
                double baseEst=-1.0*BIT_NUM2*log(((double)BIT_NUM2-lastBitmap.count())/BIT_NUM2);
                double factor=pow(2,baseIdx);
                return round(baseEst*factor);
            }
        }
        double baseEst=0;
        for(uint32_t idx=baseIdx;idx<BITMAP_NUM-1;idx++){
            baseEst-=BIT_NUM1*log(((double)BIT_NUM1-baseBitmaps[idx].count())/BIT_NUM1);
        }
        baseEst-=BIT_NUM2*log(((double)BIT_NUM2-lastBitmap.count())/BIT_NUM2);
        double factor=pow(2,baseIdx);
        return round(baseEst*factor);
    }
    void intersectOp(const MultiResolutionBitmap& MRB){
        for(uint32_t idx=0;idx<BITMAP_NUM-1;idx++){
            for(uint32_t bitIdx=0;bitIdx<BIT_NUM1;bitIdx++){
                if(!MRB.baseBitmaps[idx].test(bitIdx)){
                    baseBitmaps[idx].reset(bitIdx);
                }
            }
        }
        for(uint32_t bitIdx=0;bitIdx<BIT_NUM2;bitIdx++){
            if(!MRB.lastBitmap.test(bitIdx)){
                lastBitmap.reset(bitIdx);
            }
        }
    }
    void unionOp(const MultiResolutionBitmap& MRB){
        for(uint32_t idx=0;idx<BITMAP_NUM-1;idx++){
            for(uint32_t bitIdx=0;bitIdx<BIT_NUM1;bitIdx++){
                if(MRB.baseBitmaps[idx].test(bitIdx)){
                    baseBitmaps[idx].set(bitIdx);
                }
            }
        }
        for(uint32_t bitIdx=0;bitIdx<BIT_NUM2;bitIdx++){
            if(MRB.lastBitmap.test(bitIdx)){
                lastBitmap.set(bitIdx);
            }
        }
    }
//    uint32_t getOneBitsNum(){
//        uint32_t count=0;
//        for(uint32_t idx=0;idx<BITMAP_NUM-1;idx++){
//            for(uint32_t bitIdx=0;bitIdx<BIT_NUM1;bitIdx++){
//                if(baseBitmaps[idx].test(bitIdx)){
//                    count++;
//                }
//            }
//        }
//        for(uint32_t bitIdx=0;bitIdx<BIT_NUM2;bitIdx++){
//            if(lastBitmap.test(bitIdx)){
//                count++;
//            }
//        }
//        return count;
//    }
};

template<uint32_t BIT_NUM1,uint32_t BIT_NUM2,uint32_t SET_MAX,uint32_t SET_LAST_MAX,uint32_t BITMAP_NUM>
struct BUCKET{
    MultiResolutionBitmap<BIT_NUM1,BIT_NUM2,SET_MAX,SET_LAST_MAX,BITMAP_NUM> MRB;
    uint32_t key;
    uint8_t level;
};

template<uint32_t BIT_NUM1,uint32_t BIT_NUM2,uint32_t SET_MAX,uint32_t SET_LAST_MAX,uint32_t BITMAP_NUM,uint32_t ROW_NUM,uint32_t COL_NUM>
class SKETCH{
private:
    vector<vector<BUCKET<BIT_NUM1,BIT_NUM2,SET_MAX,SET_LAST_MAX,BITMAP_NUM>>> buckets;

public:
    SKETCH():buckets(ROW_NUM,vector<BUCKET<BIT_NUM1,BIT_NUM2,SET_MAX,SET_LAST_MAX,BITMAP_NUM>>(COL_NUM)){}
    void insert(const pair<uint32_t,uint32_t>& pkt,uint32_t HASH_SEED) {
        uint32_t hashValue1 = 0;
        MurmurHash3_x86_32(&pkt, 8, HASH_SEED, &hashValue1);
        uint32_t mszNum = __builtin_clz(hashValue1);

        for (uint32_t rowIdx = 0; rowIdx < ROW_NUM; rowIdx++) {
            uint32_t hashValue = 0;
            MurmurHash3_x86_32(&pkt.first, 4, HASH_SEEDS[rowIdx], &hashValue);
            uint32_t colIdx = hashValue % COL_NUM;
//            buckets[rowIdx][colIdx].MRB.update(pkt.first^pkt.second);
            buckets[rowIdx][colIdx].MRB.update(pkt,mszNum,HASH_SEED);

            if (buckets[rowIdx][colIdx].level <= mszNum) {
                buckets[rowIdx][colIdx].key = pkt.first;
                buckets[rowIdx][colIdx].level = mszNum;
            }
        }
    }
    uint32_t query(uint32_t key)const{
        uint32_t hashValue2 = 0;
        MurmurHash3_x86_32(&key,4,HASH_SEEDS[0], &hashValue2);
        uint32_t colIdx = hashValue2 % COL_NUM;

        MultiResolutionBitmap<BIT_NUM1,BIT_NUM2,SET_MAX,SET_LAST_MAX,BITMAP_NUM> newMRB=buckets[0][colIdx].MRB;

        for (uint32_t rowIdx = 1; rowIdx < ROW_NUM; rowIdx++) {
            MurmurHash3_x86_32(&key,4,HASH_SEEDS[rowIdx], &hashValue2);
            colIdx = hashValue2 % COL_NUM;
            newMRB.intersectOp(buckets[rowIdx][colIdx].MRB);
        }
        return newMRB.query();
    }

    void getEstimatedFlowSpreads(unordered_map<uint32_t, uint32_t>& estimatedFlowSpreads) const{
        for (uint32_t rowIdx = 0; rowIdx < ROW_NUM; rowIdx++) {
            for (uint32_t colIdx = 0; colIdx < COL_NUM; colIdx++) {
                uint32_t key=buckets[rowIdx][colIdx].key;

                auto iter2=estimatedFlowSpreads.find(key);
                if(iter2==estimatedFlowSpreads.end()){
                    estimatedFlowSpreads[key]= query(key);
                }
            }
        }
    }
//    uint32_t getOneBitsNum(){
//
//        uint32_t count=0;
//
//        for (uint32_t rowIdx = 0; rowIdx < ROW_NUM; rowIdx++) {
//            for (uint32_t colIdx = 0; colIdx < COL_NUM; colIdx++) {
//                count+=buckets[rowIdx][colIdx].MRB.getOneBitsNum();
//            }
//        }
//        return count;
//    }
    static void getMergedEstimatedFlowSpreads(const vector<SKETCH>& sketches,unordered_map<uint32_t, uint32_t>& estimatedFlowSpreads){
        SKETCH<BIT_NUM1,BIT_NUM2,SET_MAX,SET_LAST_MAX,BITMAP_NUM,ROW_NUM,COL_NUM> mergedSketch;

        for(uint32_t sketchIdx=0;sketchIdx<sketches.size();sketchIdx++){
            for (uint32_t rowIdx = 0; rowIdx < ROW_NUM; rowIdx++) {
                for (uint32_t colIdx = 0; colIdx < COL_NUM; colIdx++) {
                    mergedSketch.buckets[rowIdx][colIdx].MRB.unionOp(sketches[sketchIdx].buckets[rowIdx][colIdx].MRB);
                    if(mergedSketch.buckets[rowIdx][colIdx].level<sketches[sketchIdx].buckets[rowIdx][colIdx].level || sketchIdx==0){
                        mergedSketch.buckets[rowIdx][colIdx].level=sketches[sketchIdx].buckets[rowIdx][colIdx].level;
                        mergedSketch.buckets[rowIdx][colIdx].key=sketches[sketchIdx].buckets[rowIdx][colIdx].key;
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

template<uint32_t BIT_NUM1,uint32_t BIT_NUM2,uint32_t SET_MAX,uint32_t SET_LAST_MAX,uint32_t BITMAP_NUM,uint32_t ROW_NUM,uint32_t COL_NUM>
void processMultipleTracesAndGetAvgMetrics(const char* traceFileDir,uint32_t fileNum,uint32_t step,vector<vector<pair<uint32_t,uint32_t>>>& pktsOfFiles,vector<unordered_map<uint32_t, uint32_t>>& actualFlowSpreadsOfFiles,vector<uint32_t>& actualItemNumOfFiles,const vector<uint32_t>& topKs, string outputDirPath,string traceInfo){
    uint32_t totalMem=((BIT_NUM1*(BITMAP_NUM-1)+BIT_NUM2)+32+5)*ROW_NUM*COL_NUM;
    cout << "totalMem:" << totalMem/(8*1024.0) << "KB" << endl;
    cout << "*********************" << endl;

    char temp[500]{0};
    sprintf(temp,"%s data %d files src SS b=%d bL=%d sM=%d sML=%d c=%d r=%d w=%d mem=%.2fKB",traceInfo.c_str(),fileNum,BIT_NUM1,BIT_NUM2,SET_MAX,SET_LAST_MAX,BITMAP_NUM,ROW_NUM,COL_NUM,totalMem/(8.0*1024));
    string savedFileName(temp);

    double avgThroughput=0;
    vector<double> THs(topKs.size()),TPs(topKs.size()),FPs(topKs.size()),FNs(topKs.size()),PRs(topKs.size());
    vector<double> RCs(topKs.size()),F1s(topKs.size()),AREs(topKs.size()),AAEs(topKs.size()),realFlowNums(topKs.size()),fakeFlowNums(topKs.size());
    vector<vector<double>> avgMetrics={THs,PRs,RCs,F1s,AREs,AAEs,TPs,FPs,FNs,realFlowNums,fakeFlowNums};
    string metricsFilePath = outputDirPath+savedFileName + "-avg.metrics";
    ofstream metricsFile;
    unordered_map<double, vector<double>> avgF1;
    for(int eachTime=0;eachTime<10;eachTime++){
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
            SKETCH<BIT_NUM1,BIT_NUM2,SET_MAX,SET_LAST_MAX,BITMAP_NUM,ROW_NUM,COL_NUM> sketch;
        
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

    const uint32_t BIT_NUM1=64;
    const uint32_t BIT_NUM2=176;
    const uint32_t SET_MAX=60;
    const uint32_t SET_LAST_MAX=162;
    const uint32_t ROW_NUM=4;

    // params for 2019 single measurement point
    // 100KB
    const uint32_t BITMAP_NUM_2019_100KB=12;
    const uint32_t COL_NUM_2019_100KB=223;
    // 200KB
    const uint32_t BITMAP_NUM_2019_200KB=12;
    const uint32_t COL_NUM_2019_200KB=447;
    // 300KB
    const uint32_t BITMAP_NUM_2019_300KB=12;
    const uint32_t COL_NUM_2019_300KB=670;
    // 150KB
    const uint32_t BITMAP_NUM_2019_150KB=12;
    const uint32_t COL_NUM_2019_150KB=335;
    // 250KB
    const uint32_t BITMAP_NUM_2019_250KB=12;
    const uint32_t COL_NUM_2019_250KB=558;

    const char* traceFileDir2019=R"(../data/2019ipv4/)";
    
    uint32_t fileNum=5;
    uint32_t step=5;
    vector<vector<pair<uint32_t,uint32_t>>> pktsOfFiles;
    vector<unordered_map<uint32_t, uint32_t>> actualFlowSpreadsOfFiles;
    vector<uint32_t> actualItemNumOfFiles;
    getDataSets(traceFileDir2019,fileNum,step,pktsOfFiles,actualFlowSpreadsOfFiles,actualItemNumOfFiles);
    cout<<"get data set finished"<<endl;
    
    
//     //single measurement point
    processMultipleTracesAndGetAvgMetrics<BIT_NUM1,BIT_NUM2,SET_MAX,SET_LAST_MAX,BITMAP_NUM_2019_100KB,ROW_NUM,COL_NUM_2019_100KB>(traceFileDir2019,fileNum,step,pktsOfFiles,actualFlowSpreadsOfFiles,actualItemNumOfFiles,topKs,outputDirPath,"2019");
    processMultipleTracesAndGetAvgMetrics<BIT_NUM1,BIT_NUM2,SET_MAX,SET_LAST_MAX,BITMAP_NUM_2019_200KB,ROW_NUM,COL_NUM_2019_200KB>(traceFileDir2019,fileNum,step,pktsOfFiles,actualFlowSpreadsOfFiles,actualItemNumOfFiles,topKs,outputDirPath,"2019");
    processMultipleTracesAndGetAvgMetrics<BIT_NUM1,BIT_NUM2,SET_MAX,SET_LAST_MAX,BITMAP_NUM_2019_300KB,ROW_NUM,COL_NUM_2019_300KB>(traceFileDir2019,fileNum,step,pktsOfFiles,actualFlowSpreadsOfFiles,actualItemNumOfFiles,topKs,outputDirPath,"2019");
    processMultipleTracesAndGetAvgMetrics<BIT_NUM1,BIT_NUM2,SET_MAX,SET_LAST_MAX,BITMAP_NUM_2019_150KB,ROW_NUM,COL_NUM_2019_150KB>(traceFileDir2019,fileNum,step,pktsOfFiles,actualFlowSpreadsOfFiles,actualItemNumOfFiles,topKs,outputDirPath,"2019");
    processMultipleTracesAndGetAvgMetrics<BIT_NUM1,BIT_NUM2,SET_MAX,SET_LAST_MAX,BITMAP_NUM_2019_150KB,ROW_NUM,COL_NUM_2019_250KB>(traceFileDir2019,fileNum,step,pktsOfFiles,actualFlowSpreadsOfFiles,actualItemNumOfFiles,topKs,outputDirPath,"2019");

    return 0;
}
