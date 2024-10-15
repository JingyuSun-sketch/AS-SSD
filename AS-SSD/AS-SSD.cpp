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
const uint32_t HASH_SEEDS[]={2725628290,1444492164,1001369998,2291192083,3632694421,3217476757,962572696,3523028760,1564808096,1686744101};
uint32_t R[32]={2572712615,186525479,1267589045,2944996661,3807606838,2892595386,1879920448,3800713281,4044617030,1088214349,2912533710,4048616015,223249368,1094862426,1603814236,3377932768,121561825,2119485156,1642325476,1361943285,3238649846,992916862,860938402,674876105,3044024791,706198233,1838710314,2602646517,2912533710,4048616015,223249368,1094862426};

string avgF1RatioFileName = "./test/F1OfRatios.txt";//

ofstream avgF1RatioFile;

bool cmpPairFunc(pair<uint32_t,uint32_t>p1, pair<uint32_t,uint32_t>p2)
{
    return p1.second > p2.second;
}

class BucketArrayFreeRs{
    private:
        uint32_t bucketNum;
        uint32_t regNum;
        uint32_t bucketCellNum;
        double q_r=1.0;
        
        vector<vector<pair<uint32_t,uint32_t>>> bucketArray;
        uint8_t* regArray;
        
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

        int updateBucket(const pair<uint32_t,uint32_t>&pkt,uint32_t hashBucketIndex,double ap_i,double q_r,double sigmas,uint32_t HASH_SEED){
            uint32_t minCellIdx=-1,minCellVal=-1;
                 
            for(int cellIdx=0;cellIdx<bucketCellNum;cellIdx++){
                if(bucketArray[hashBucketIndex][cellIdx].second==0){
                    bucketArray[hashBucketIndex][cellIdx].first=pkt.first;
                    
                    double updateValue=1/q_r;
                    uint32_t updateValueToInt = floor(updateValue);

                    double temp=(double)rand()/RAND_MAX;
                    if(temp<(updateValue- updateValueToInt)){
                        updateValueToInt+=1;
                    }
                    
                    bucketArray[hashBucketIndex][cellIdx].second=updateValueToInt;
                    return 1;
                }
                else if(bucketArray[hashBucketIndex][cellIdx].first==pkt.first){ 
                    double hash_adaptive_index=(double)2.0*bucketArray[hashBucketIndex][cellIdx].second*sigmas+1;                   
                    double adaptive_probaility=1.0/max(ap_i,hash_adaptive_index);
                    
                    if(q_r<adaptive_probaility){
                        double updateValue=1/q_r;
                        uint32_t updateValueToInt = floor(updateValue);

                        double temp=(double)rand()/RAND_MAX;
                        if(temp<(updateValue- updateValueToInt)){
                            updateValueToInt+=1;
                        }
                        bucketArray[hashBucketIndex][cellIdx].second+=updateValueToInt;
                        return 1;
                    }else{
                        uint32_t hash_value4=0;
                        MurmurHash3_x86_32(&pkt,8,HASH_SEED+91323187,&hash_value4);
                        double adaptive_num = (double)hash_value4/(pow(2,32));
                        
                        if(adaptive_num>=adaptive_probaility/q_r){
                            return 0;
                        }
                        
                        double updateValue=1/adaptive_probaility;
                        uint32_t updateValueToInt = floor(updateValue);

                        double temp=(double)rand()/RAND_MAX;
                        if(temp<(updateValue- updateValueToInt)){
                            updateValueToInt+=1;
                        }

                        bucketArray[hashBucketIndex][cellIdx].second+=updateValueToInt;
                        return 1;
                    }
                }
                else{
                    if(bucketArray[hashBucketIndex][cellIdx].second<minCellVal){
                        minCellIdx=cellIdx;
                        minCellVal=bucketArray[hashBucketIndex][cellIdx].second;
                    }
                }
            }
            
            double hash_adaptive_index=(double)2.0*bucketArray[hashBucketIndex][minCellIdx].second*sigmas+1;                   
            double adaptive_probaility=1.0/max(ap_i,hash_adaptive_index);
            
            double p_t=min(1.0,adaptive_probaility/q_r);
            
            uint32_t hash_value5=0;
            MurmurHash3_x86_32(&pkt,8,HASH_SEED+25743187,&hash_value5);
            double temp1 = (double)hash_value5/(pow(2,32));
            
            if(temp1<p_t){
                
                double updateValue=1/q_r;
                uint32_t updateValueToInt = floor(updateValue);
                double temp=(double)rand()/RAND_MAX;
                if(temp<(updateValue- updateValueToInt)){
                    updateValueToInt+=1;
                }
                
                double changeBucketTP=1.0/((minCellVal+updateValueToInt)*(q_r*p_t));
                double temp2=(double)rand()/RAND_MAX;
                if(temp2<changeBucketTP){
                    bucketArray[hashBucketIndex][minCellIdx].first=pkt.first;
                    bucketArray[hashBucketIndex][minCellIdx].second=minCellVal+updateValueToInt;
                }
                return 1;
            }
            return 0;
        }

    public:
        BucketArrayFreeRs(uint32_t bucketNum,uint32_t regNum,uint32_t bucketCellNum):bucketArray(bucketNum,vector<pair<uint32_t,uint32_t>>(bucketCellNum)){
            this->bucketNum=bucketNum;
            this->regNum=regNum;
            this->bucketCellNum=bucketCellNum;
            regArray=new uint8_t[int(ceil(regNum*5.0/8.0))]();
        }
        void updateBucketArray(const pair<uint32_t,uint32_t>& pkt,double ap_i,uint32_t HASH_SEED){
            uint32_t key=pkt.first;
            uint32_t element=pkt.second;

            uint32_t hash_value=0;
            MurmurHash3_x86_32(&pkt,8,HASH_SEED,&hash_value);
            uint32_t hash_index=hash_value%regNum;

            uint32_t hash_value2=0;
            MurmurHash3_x86_32(&pkt,8,HASH_SEED+15416691,&hash_value2);
            uint32_t p_e=min(31,__builtin_clz(hash_value2)+1);

            uint32_t hash_value3=0;
            MurmurHash3_x86_32(&pkt.first,4,HASH_SEED+15417891,&hash_value3);
            uint32_t hashBucketIndex=hash_value3%bucketNum;
            // cout<<"hash_value3:"<<hash_value3<<endl;

            // cout<<"emptyIndex:"<<emptyIndex<<endl;
            uint32_t regValue=getRegValue(hash_index);
            // cout<<"regValue:"<<regValue<<endl;
            if(p_e > regValue){
                double sigma=0.05;
                double sigmas=pow(sigma,2);
                int charge=updateBucket(pkt,hashBucketIndex,ap_i,q_r,sigmas,HASH_SEED);
                if(charge==0){
                    return;
                }
                if(p_e != 31){
                    q_r+=(pow(2,-1.0*p_e)-pow(2,-1.0*regValue))/regNum;
                }else{
                    q_r-=pow(2,-1.0*regValue)/regNum;
                }
                setRegValue(hash_index,p_e);
                
                
                // cout<<"set"<<endl;
            }
        }

        void getEstimatedFlowSpreads(unordered_map<uint32_t, uint32_t>& estimatedFlowSpreads)const{
            for (uint32_t bucketidx=0;bucketidx<bucketNum;bucketidx++){
                for(int cellIdx=0;cellIdx<bucketCellNum;cellIdx++){
                    uint32_t key=bucketArray[bucketidx][cellIdx].first;
                    uint32_t estSpread=bucketArray[bucketidx][cellIdx].second;

                    estimatedFlowSpreads[key]=estSpread;
                }
            }
        }
    
        double getQR(){
            return q_r;
        }
        
        ~BucketArrayFreeRs(){
            delete []regArray;
        }
};

// bool cmpPairFunc(pair<uint32_t,uint32_t>p1, pair<uint32_t,uint32_t>p2)
// {
//     return p1.second > p2.second;
// }

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

void processMultipleTracesAndGetAvgMetrics(double ap_j,uint32_t bucketNum,uint32_t bucketCellNum,uint32_t regNum,const char* traceFileDir,uint32_t fileNum,uint32_t step,vector<vector<pair<uint32_t,uint32_t>>>& pktsOfFiles,vector<unordered_map<uint32_t, uint32_t>>& actualFlowSpreadsOfFiles,vector<uint32_t>& actualItemNumOfFiles,const vector<uint32_t>& topKs, string outputDirPath,string traceInfo)
{
    //bucketNum
    //bucketCellsize
    //regNum

    uint32_t totalMem=bucketNum*bucketCellNum*64+regNum*5;
    double bucketMemRatio=(double)bucketCellNum*bucketNum*64/totalMem;
    cout << "totalMem:" << totalMem/(8*1024.0) << "KB" << endl;
    cout << "*********************" << endl;

    char temp[500]{0};
    sprintf(temp,"New %s data %d files %d step src Bkt_AP_MaxHalf_freeRS ap=%.2f sigma=%.2f rGNum=%d bktNum=%d bktSz=%d bktR=%.2f mem=%.2fKB",traceInfo.c_str(),fileNum,step,ap_j,0.05,regNum,bucketNum,bucketCellNum,bucketMemRatio,totalMem/(8.0*1024));
    string savedFileName(temp);
    double avgThroughput=0;
    vector<double> THs(topKs.size()),TPs(topKs.size()),FPs(topKs.size()),FNs(topKs.size()),PRs(topKs.size());
    vector<double> RCs(topKs.size()),F1s(topKs.size()),AREs(topKs.size()),AAEs(topKs.size()),realFlowNums(topKs.size()),fakeFlowNums(topKs.size());
    vector<vector<double>> avgMetrics={THs,PRs,RCs,F1s,AREs,AAEs,TPs,FPs,FNs,realFlowNums,fakeFlowNums};
    string metricsFilePath = outputDirPath+savedFileName + "-avg.metrics";
    ofstream metricsFile;
    unordered_map<double, vector<double>> avgF1;
    for(int eachTime=0;eachTime<10;eachTime++){
        uint32_t HASH_SEED = HASH_SEEDS[eachTime];
        cout<<HASH_SEED<<endl;
        for(uint32_t i=step;i<=fileNum;i=i+step){
            char traceFilePath[256]{0};
            sprintf(traceFilePath,"%s%02d.dat--%02d.dat",traceFileDir,i-step,i-1);
            BucketArrayFreeRs sketch(bucketNum,regNum,bucketCellNum);
            unordered_map<uint32_t, uint32_t> actualFlowSpreads=actualFlowSpreadsOfFiles[(i/step)-1];
            unsigned int actualItemNum=0;
            
            clock_t time1 = clock();
            for(uint32_t j=(i-step);j<i;j++){
                    //prepare dataset
                vector<pair<uint32_t,uint32_t>> pkts=pktsOfFiles[j];
                // unordered_map<uint32_t, uint32_t> actualFlowSpreadsPart=actualFlowSpreadsOfFiles[j];
                unsigned int actualItemNumPart = actualItemNumOfFiles[j];
                cout<<"pktsSize"<<pkts.size()<<endl;
                cout << "prepare algorithm "<< endl;
                // clock_t time1 = clock();
                for (unsigned int i = 0; i < pkts.size(); i++) {
                    sketch.updateBucketArray(pkts[i],ap_j,HASH_SEED);
                }
                
                // metricsFile.open(metricsFilePath,ios::app);
                // double qr=sketch.getQR();
                // metricsFile<<"min:"<<j<<" q_r:"<<qr<<endl;
                // metricsFile.close();
                
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

    //         string metricsFilePath = outputDirPath+savedFileName + "-avg.metrics";
    //         ofstream metricsFile;
            metricsFile.open(metricsFilePath,ios::app);
            metricsFile<<"throughput: " << throughput << " Mpps" << ", each insert operation uses " << 1000.0 / throughput << " ns" << endl;
            metricsFile<<"ap:"<<ap_j<<endl;
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

    vector<double> ap_js={1.0};
    string outputDirPath="./test/results2019MemTen/";//test
    const char* traceFileDir=R"(../data/2019ipv4/)";
//     const char* traceFileDir=R"(../data/Mul/)";
//     const char* traceFileDir=R"(/home/tmp/Desktop/ZBY/throughput/perSource_LittleEndian/data/Mul2020.4/)";
//     const char* traceFileDir=R"(/home/tmp/Desktop/ZBY/throughput/perSource_LittleEndian/data/Cos/)";

    uint32_t fileNum=10;
    uint32_t step=10;
    vector<vector<pair<uint32_t,uint32_t>>> pktsOfFiles;
    vector<unordered_map<uint32_t, uint32_t>> actualFlowSpreadsOfFiles;
    vector<uint32_t> actualItemNumOfFiles;

    getDataSets(traceFileDir,fileNum,step,pktsOfFiles,actualFlowSpreadsOfFiles,actualItemNumOfFiles);
    cout<<"get data set finished"<<endl;
    uint32_t bucketCellSize=8;
    for(double ap_j:ap_js){
        for(int mem=100;mem<=150;mem+=50){

            avgF1RatioFile.open(avgF1RatioFileName,ios::app);
            avgF1RatioFile<<"mem="<<mem<<" KB"<<endl;
            avgF1RatioFile.close();

            for(int i=45;i<=80;i+=5){
                double BK_FreeRSMemRatio=i*0.01;
                uint32_t totalMemBits=mem*8*1024;
                double BK_FreeRSMem=totalMemBits*BK_FreeRSMemRatio;
                uint32_t BK_FreeRSSize=uint32_t(BK_FreeRSMem/64.0);
                int bucketNum = floor(BK_FreeRSSize/bucketCellSize);
                int cellMore = BK_FreeRSSize%bucketCellSize;
                int halfBucketCellSize=bucketCellSize/2;
                if(cellMore>halfBucketCellSize){
                    bucketNum+=1;
                }
                
                double regMem=totalMemBits-bucketNum*bucketCellSize*64;
                uint32_t regNum=uint32_t(regMem/5.0);

                processMultipleTracesAndGetAvgMetrics(ap_j,bucketNum,bucketCellSize,regNum,traceFileDir,fileNum,step,pktsOfFiles,actualFlowSpreadsOfFiles,actualItemNumOfFiles,topKs,outputDirPath,"2019");

            }
            avgF1RatioFile.open(avgF1RatioFileName,ios::app);
            avgF1RatioFile<<endl;
            avgF1RatioFile.close();
        }
    }
    return 0;  
}
