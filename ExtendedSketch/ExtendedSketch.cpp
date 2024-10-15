#include <cstdint>
#include <cstring>
#include <iostream>
#include <cmath>
#include <random>
#include <cstdlib>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <ctime>
#include <vector>
#include <iomanip>
#include <numeric>
#include "MurmurHash3.h"
#include <fstream>
#include <string>
#include "mini-gmp.h"

using namespace std;
const uint32_t HASH_SEEDS[]={2725628290,1444492164,1001369998,2291192083,3632694421,3217476757,962572696,3523028760,1564808096,1686744101};

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

template<uint32_t BASE_BIT_NUM>
class BUCKET{
public:
    BITMAP EC;
    unordered_set<uint32_t> flag;
    uint8_t ET=0;
    BUCKET():EC(BASE_BIT_NUM){};
    void update(const pair<uint32_t,uint32_t>& pkt,uint32_t HASH_SEED){
        uint32_t hashValue1=0;
        MurmurHash3_x86_32(&pkt,8,HASH_SEED, &hashValue1);

        uint32_t bitIdx=hashValue1%(BASE_BIT_NUM<<ET);
        EC.set(bitIdx);
    }
    void update(const pair<uint32_t,uint32_t>& pkt,uint32_t nextColIdx,uint32_t HASH_SEED){
        uint32_t hashValue1=0;
        MurmurHash3_x86_32(&pkt,8,HASH_SEED, &hashValue1);

        uint32_t bitIdx=hashValue1%(BASE_BIT_NUM<<ET);
        EC.set(bitIdx);

//        flag.insert(nextColIdx);
        if(flag.find(nextColIdx)==flag.end()){
            flag.insert(nextColIdx);
        }

    }
    void extend(uint32_t rowIdx){
        uint32_t curBitNum=BASE_BIT_NUM<<ET;
        uint32_t newBitNum=BASE_BIT_NUM<<(ET+1);

        BITMAP newEC(newBitNum);

        if(rowIdx%2==0){
            for(uint32_t bitIdx=0;bitIdx<curBitNum;bitIdx++){
                if(EC.test(bitIdx)){
                    newEC.set(curBitNum+bitIdx);
                }
            }
        }else{
            for(uint32_t bitIdx=0;bitIdx<curBitNum;bitIdx++){
                if(EC.test(bitIdx)){
                    newEC.set(bitIdx);
                }
            }
        }
        EC=newEC;
        ET++;
    }
    double getOneBitsRatio() const{
        return EC.getOneBitsRatio();
    }
    uint32_t query(double coe) const{
        uint32_t curBitNum=BASE_BIT_NUM<<ET;
        double emptyRatio=EC.getZeroBitsRatio();

        if(emptyRatio==0){
            emptyRatio=1.0/curBitNum;
        }
        double baseVal=-1.0*curBitNum*log(emptyRatio);

        uint32_t sum=0;
        for(uint32_t i=0;i<ET;i++){
            sum+=1<<i;
        }
        double term=coe*BASE_BIT_NUM*sum;

        return round(baseVal+term);
    }
};

template<uint32_t ROW_NUM, uint32_t BASE_BIT_NUM, uint8_t ET_TH>
class SKETCH{
private:
    vector<vector<BUCKET<BASE_BIT_NUM>>> buckets;
    vector<uint32_t> COL_NUMS;
    double epsilon;

    //parameters for recover keys
    double coe;
    mpz_t P;
    mpz_t Q[ROW_NUM];
    mpz_t Q_inv[ROW_NUM];

public:
    SKETCH(vector<uint32_t>& colnums,double epsilon):COL_NUMS(colnums),buckets(ROW_NUM){
        this->epsilon=epsilon;
        for(uint32_t i=0;i<ROW_NUM;i++){
            buckets[i].resize(COL_NUMS[i]);
        }

        coe=log(pow(2.0-epsilon,2)/(4.0*(1-epsilon)));

        mpz_init_set_ui(P,1);
        for(uint32_t i=0;i<ROW_NUM;i++){
            mpz_mul_ui(P,P,COL_NUMS[i]);
        }
        for(uint32_t i=0;i<ROW_NUM;i++){
            mpz_init(Q[i]);
            mpz_tdiv_q_ui(Q[i],P,COL_NUMS[i]);
            mpz_init(Q_inv[i]);

            mpz_t colNum;
            mpz_init_set_ui(colNum,COL_NUMS[i]);
            getModularInverses(Q_inv[i],Q[i],colNum);
            mpz_clear(colNum);
        }
    }

    SKETCH(const SKETCH& other):COL_NUMS(other.COL_NUMS),buckets(other.buckets){
        epsilon=other.epsilon;
        coe=other.coe;

        mpz_init_set_ui(P,1);
        for(uint32_t i=0;i<ROW_NUM;i++){
            mpz_mul_ui(P,P,COL_NUMS[i]);
        }
        for(uint32_t i=0;i<ROW_NUM;i++){
            mpz_init(Q[i]);
            mpz_tdiv_q_ui(Q[i],P,COL_NUMS[i]);
            mpz_init(Q_inv[i]);

            mpz_t colNum;
            mpz_init_set_ui(colNum,COL_NUMS[i]);
            getModularInverses(Q_inv[i],Q[i],colNum);
            mpz_clear(colNum);
        }
    }

    ~SKETCH(){
        mpz_clear(P);

        for(int i=0;i<ROW_NUM;i++){
            mpz_clear(Q[i]);
            mpz_clear(Q_inv[i]);
        }
    }

    void insert(const pair<uint32_t,uint32_t>& pkt,uint32_t HASH_SEED) {
        uint32_t colIdxes[ROW_NUM]{0};
        uint32_t fullBucketNum=0;
        for (uint32_t rowIdx = 0; rowIdx < ROW_NUM; rowIdx++) {
            uint32_t colIdx=pkt.first%COL_NUMS[rowIdx];
            colIdxes[rowIdx]=colIdx;

            if(buckets[rowIdx][colIdx].getOneBitsRatio()>=epsilon){
                fullBucketNum++;
            }
        }

        for (uint32_t rowIdx = 0; rowIdx < ROW_NUM; rowIdx++) {
            uint32_t colIdx=colIdxes[rowIdx];
            BUCKET<BASE_BIT_NUM>& bucket=buckets[rowIdx][colIdx];

            if(fullBucketNum>=ROW_NUM){
                bucket.extend(rowIdx);
            }

            if(rowIdx!=ROW_NUM-1){
                bucket.update(pkt,colIdxes[rowIdx+1],HASH_SEED);
            }else{
                bucket.update(pkt,HASH_SEED);
            }
        }
    }

    uint32_t query(uint32_t key) const{
        bool isExtended=false;
        for (uint32_t rowIdx = 0; rowIdx < ROW_NUM; rowIdx++) {
            uint32_t colIdx = key % COL_NUMS[rowIdx];
            if(buckets[rowIdx][colIdx].ET>0){
                isExtended=true;
                break;
            }
        }

        if(isExtended){
            return query_extension(key);
        }else{
            return query_no_extension(key);
        }
    }

    bool getEstimatedFlowSpreads(unordered_map<uint32_t, uint32_t>& estimatedFlowSpreads){
        vector<uint32_t> idxList;
        uint32_t flowCount=0;
        for(uint32_t colIdx=0;colIdx<COL_NUMS[0];colIdx++){
            if(getKeyValueFromBucket(0,colIdx,idxList,estimatedFlowSpreads,flowCount)==false){
                return false;
            }
        }
        return true;
    }

    static bool getMergedEstimatedFlowSpreads(const vector<SKETCH>& sketches,unordered_map<uint32_t, uint32_t>& estimatedFlowSpreads){
        vector<uint32_t> COL_NUMS=sketches[0].COL_NUMS;
        double epsilon=sketches[0].epsilon;
        SKETCH<ROW_NUM,BASE_BIT_NUM,ET_TH> mergedSketch(COL_NUMS,epsilon);

        for (uint32_t rowIdx = 0; rowIdx < ROW_NUM; rowIdx++) {
            for (uint32_t colIdx = 0; colIdx < COL_NUMS[rowIdx]; colIdx++) {
                BUCKET<BASE_BIT_NUM>& mergedBucket=mergedSketch.buckets[rowIdx][colIdx];
                for(uint32_t sketchIdx=0;sketchIdx<sketches.size();sketchIdx++){
                    BUCKET<BASE_BIT_NUM> bucket=sketches[sketchIdx].buckets[rowIdx][colIdx];

                    while(mergedBucket.ET<bucket.ET){
                        mergedBucket.extend(rowIdx);
                    }
                    while(bucket.ET<mergedBucket.ET){
                        bucket.extend(rowIdx);
                    }
                    uint32_t curBitNum=BASE_BIT_NUM<<bucket.ET;
                    for(uint32_t bitIdx=0;bitIdx<curBitNum;bitIdx++){
                        if(bucket.EC.test(bitIdx)){
                            mergedBucket.EC.set(bitIdx);
                        }
                    }

                    for(auto iter=bucket.flag.begin();iter!=bucket.flag.end();iter++){
                        mergedBucket.flag.insert(*iter);
                    }
                }

            }
        }

        if(mergedSketch.getEstimatedFlowSpreads(estimatedFlowSpreads)==false){
            return false;
        }else{
            return true;
        }
    }
private:

    uint32_t query_no_extension(uint32_t key) const{
        uint32_t zeroNum=0;

        uint32_t colIdxes[ROW_NUM]{0};
        for (uint32_t rowIdx = 0; rowIdx < ROW_NUM; rowIdx++) {
            uint32_t colIdx=key%COL_NUMS[rowIdx];
            colIdxes[rowIdx]=colIdx;
        }

        for(uint32_t bitIdx=0;bitIdx<BASE_BIT_NUM;bitIdx++){
            for (uint32_t rowIdx = 0; rowIdx < ROW_NUM; rowIdx++) {
                if(buckets[rowIdx][colIdxes[rowIdx]].EC.test(bitIdx)==false){
                    zeroNum++;
                    break;
                }
            }
        }

        if(zeroNum!=0){
            double val=BASE_BIT_NUM*log((double)BASE_BIT_NUM/zeroNum);
            return round(val);
        }else{
            double val=BASE_BIT_NUM*log(BASE_BIT_NUM);
            return round(val);
        }
    }

    uint32_t query_extension(uint32_t key) const{
        uint32_t minVal=-1;
        for (uint32_t rowIdx = 0; rowIdx < ROW_NUM; rowIdx++) {
            uint32_t colIdx=key%COL_NUMS[rowIdx];
            uint32_t val=buckets[rowIdx][colIdx].query(coe);
            if(val<minVal){
                minVal=val;
            }
        }
        return minVal;
    }

    void getModularInverses(mpz_t& result,mpz_t a,mpz_t m){
        mpz_t x,y,g;
        mpz_init(x);
        mpz_init(y);
        mpz_init(g);

        mpz_gcdext(g,x,y,a,m);

        if(mpz_cmp_si(g,1)!=0){
            mpz_set_ui(result,0);
        }else{
            if(mpz_cmp_si(x,0)<0){
                mpz_tdiv_r(x,x,m);
                mpz_add(x,x,m);
                mpz_tdiv_r(x,x,m);
            }
            mpz_set(result,x);
        }
        mpz_clear(x);
        mpz_clear(y);
        mpz_clear(g);
    }

    uint32_t reversibleCalculation(vector<uint32_t>& idxList,bool& result) const{
        mpz_t sum,temp;
        mpz_init_set_ui(sum,0);
        for(uint32_t rowIdx=0;rowIdx<ROW_NUM;rowIdx++){
            mpz_init_set_ui(temp,1);
            mpz_mul(temp,temp,Q[rowIdx]);
            mpz_mul(temp,temp,Q_inv[rowIdx]);
            mpz_mul_ui(temp,temp,idxList[rowIdx]);
            mpz_add(sum,sum,temp);
            mpz_clear(temp);
        }

        mpz_init_set_ui(temp,0);
        mpz_tdiv_r(temp,sum,P);

        int cmpResult=mpz_cmp_ui(temp,UINT32_MAX);
        if(cmpResult<=0){
            result=true;
        }else{
            result=false;
        }

        mpz_t mask;
        mpz_init_set_ui(mask,UINT32_MAX);
        mpz_and(temp,temp,mask);

        uint32_t key= mpz_get_ui(temp);

        mpz_clear(sum);
        mpz_clear(mask);
        mpz_clear(temp);

        return key;
    }


    //不保存key，得到一个直接估计
    bool getKeyValueFromBucket(uint32_t rowIdx,uint32_t colIdx,vector<uint32_t>& idxList,unordered_map<uint32_t, uint32_t>& estimatedFlowSpreads,uint32_t& flowCount){
//         if(flowCount>0 && flowCount%1000000==0){
//             cout<<"get "<<flowCount<<" keys, record "<<estimatedFlowSpreads.size()<<"keys"<<endl;
//         }
        if(flowCount>1000000000){//too many fake flows
            cout<<"too many fake flows "<<endl;
            return false;
        }
        const BUCKET<BASE_BIT_NUM>& bucket=buckets[rowIdx][colIdx];

        if(bucket.ET>=ET_TH){
            idxList.push_back(colIdx);
            if(rowIdx>=ROW_NUM-1){
                bool isGetKey=true;
                uint32_t key= reversibleCalculation(idxList,isGetKey);
                flowCount++;

                if(isGetKey){
                    uint32_t estSpread=query(key);
                    estimatedFlowSpreads[key]=estSpread;
                }
            }else{
                const unordered_set<uint32_t>& flag=bucket.flag;
                for(auto iter=flag.begin();iter!=flag.end();iter++){
                    uint32_t nextColIdx=*iter;
                    if(getKeyValueFromBucket(rowIdx+1,nextColIdx,idxList,estimatedFlowSpreads,flowCount)==false){
                        return false;
                    }
                }
            }
            idxList.pop_back();
        }
        return true;
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

    for(uint32_t k:topKs){
        double threshold=0.5*(actualFlowSpreadsVec[k-1].second+actualFlowSpreadsVec[k].second);

        unordered_set<uint32_t> trueSuperSpreads;
        for(uint32_t i=0;i<actualFlowSpreadsVec.size();i++){
            if(actualFlowSpreadsVec[i].second>=threshold){
                trueSuperSpreads.insert(actualFlowSpreadsVec[i].first);
            }
        }

        double TP=0,FP=0,totalRE=0,totalAE=0,fakeFlowNum=0;

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

template<uint32_t ROW_NUM, uint32_t BASE_BIT_NUM, uint8_t ET_TH>
void processMultipleTracesAndGetAvgMetrics(const char* traceFileDir,uint32_t fileNum,uint32_t step,vector<vector<pair<uint32_t,uint32_t>>>& pktsOfFiles,vector<unordered_map<uint32_t, uint32_t>>& actualFlowSpreadsOfFiles,vector<uint32_t>& actualItemNumOfFiles,const vector<uint32_t>& topKs, string outputDirPath,string traceInfo,vector<uint32_t>& COL_NUMS,double epsilon){
    uint32_t totalBucketNum=0;
    for(uint32_t colNum:COL_NUMS){
        totalBucketNum+=colNum;
    }
    //the extensible bitmap in EC has baseBitsNum bits initially, recording the number of used bits in bitmap and the number of extended times use 32 bits and 8 bits separately
    //the pointer pointing to the extensible bitmap uses 32bits, and the pointer pointing to the flag list(set) also uses 32 bits
    //the memory of the flag list is not calculated here
    uint32_t totalMem=totalBucketNum*(BASE_BIT_NUM+32+8+32+32);
    cout << "totalMem:" << totalMem/(8*1024.0) << "KB" << endl;
    cout << "*********************" << endl;

    char temp[500]{0};
    sprintf(temp,"%s data %d files src ES r=%d b=%d ep=%.2f etTh=%d mem=%.2fKB",traceInfo.c_str(),fileNum,ROW_NUM,BASE_BIT_NUM,epsilon,ET_TH,totalMem/(8.0*1024));
    string savedFileName(temp);

    double avgThroughput=0;
    vector<double> THs(topKs.size()),TPs(topKs.size()),FPs(topKs.size()),FNs(topKs.size()),PRs(topKs.size());
    vector<double> RCs(topKs.size()),F1s(topKs.size()),AREs(topKs.size()),AAEs(topKs.size()),realFlowNums(topKs.size()),fakeFlowNums(topKs.size());
    vector<vector<double>> avgMetrics={THs,PRs,RCs,F1s,AREs,AAEs,TPs,FPs,FNs,realFlowNums,fakeFlowNums};
    string metricsFilePath = outputDirPath+savedFileName + "-avg.metrics";
    ofstream metricsFile;
    unordered_map<double, vector<double>> avgF1;
    for(int eachTime=0;eachTime<1;eachTime++){
        uint32_t HASH_SEED = HASH_SEEDS[eachTime];
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
            SKETCH<ROW_NUM,BASE_BIT_NUM,ET_TH> sketch(COL_NUMS,epsilon);

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
            if(sketch.getEstimatedFlowSpreads(estFlowSpreads)==false){
                return;
            }

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
    const uint32_t BASE_BIT_NUM=32;
    const uint32_t ROW_NUM=4;

    //250KB
    vector<uint32_t> COL_NUMS_250KB{3769,3761,3767,3779};
    //150KB
    vector<uint32_t> COL_NUMS_150KB{2243,2251,2267,2273};
    //300KB
    vector<uint32_t> COL_NUMS_300KB{4517, 4513, 4547, 4493};
    //200KB
    vector<uint32_t> COL_NUMS_200KB{3011, 3001, 2999, 3037};
    //100KB
    vector<uint32_t> COL_NUMS_100KB{1499, 1511, 1489, 1523};

    uint32_t fileNum=5;
    uint32_t step=5;
    const char* traceFileDir2019=R"(../data/2019ipv4/)";
    vector<vector<pair<uint32_t,uint32_t>>> pktsOfFiles;
    vector<unordered_map<uint32_t, uint32_t>> actualFlowSpreadsOfFiles;
    vector<uint32_t> actualItemNumOfFiles;
    getDataSets(traceFileDir2019,fileNum,step,pktsOfFiles,actualFlowSpreadsOfFiles,actualItemNumOfFiles);
    cout<<"get data set finished"<<endl;

    //single measurement point
    //10 min
//    processMultipleTracesAndGetAvgMetrics<ROW_NUM,BASE_BIT_NUM,7>(traceFileDir2019,fileNum,step,pktsOfFiles,actualFlowSpreadsOfFiles,actualItemNumOfFiles,topKs,outputDirPath,"2019",COL_NUMS_300KB,epsilon);
//    processMultipleTracesAndGetAvgMetrics<ROW_NUM,BASE_BIT_NUM,7>(traceFileDir2019,fileNum,step,pktsOfFiles,actualFlowSpreadsOfFiles,actualItemNumOfFiles,topKs,outputDirPath,"2019",COL_NUMS_250KB,epsilon);
//    processMultipleTracesAndGetAvgMetrics<ROW_NUM,BASE_BIT_NUM,8>(traceFileDir2019,fileNum,step,pktsOfFiles,actualFlowSpreadsOfFiles,actualItemNumOfFiles,topKs,outputDirPath,"2019",COL_NUMS_200KB,epsilon);
//    processMultipleTracesAndGetAvgMetrics<ROW_NUM,BASE_BIT_NUM,9>(traceFileDir2019,fileNum,step,pktsOfFiles,actualFlowSpreadsOfFiles,actualItemNumOfFiles,topKs,outputDirPath,"2019",COL_NUMS_150KB,epsilon);
//    processMultipleTracesAndGetAvgMetrics<ROW_NUM,BASE_BIT_NUM,10>(traceFileDir2019,fileNum,step,pktsOfFiles,actualFlowSpreadsOfFiles,actualItemNumOfFiles,topKs,outputDirPath,"2019",COL_NUMS_100KB,epsilon);

   //1min
//    processMultipleTracesAndGetAvgMetrics<ROW_NUM,BASE_BIT_NUM,4>(traceFileDir2019,fileNum,step,pktsOfFiles,actualFlowSpreadsOfFiles,actualItemNumOfFiles,topKs,outputDirPath,"2019",COL_NUMS_300KB,epsilon);
//    processMultipleTracesAndGetAvgMetrics<ROW_NUM,BASE_BIT_NUM,4>(traceFileDir2019,fileNum,step,pktsOfFiles,actualFlowSpreadsOfFiles,actualItemNumOfFiles,topKs,outputDirPath,"2019",COL_NUMS_250KB,epsilon);
//    processMultipleTracesAndGetAvgMetrics<ROW_NUM,BASE_BIT_NUM,5>(traceFileDir2019,fileNum,step,pktsOfFiles,actualFlowSpreadsOfFiles,actualItemNumOfFiles,topKs,outputDirPath,"2019",COL_NUMS_200KB,epsilon);
//    processMultipleTracesAndGetAvgMetrics<ROW_NUM,BASE_BIT_NUM,5>(traceFileDir2019,fileNum,step,pktsOfFiles,actualFlowSpreadsOfFiles,actualItemNumOfFiles,topKs,outputDirPath,"2019",COL_NUMS_150KB,epsilon);
//    processMultipleTracesAndGetAvgMetrics<ROW_NUM,BASE_BIT_NUM,6>(traceFileDir2019,fileNum,step,pktsOfFiles,actualFlowSpreadsOfFiles,actualItemNumOfFiles,topKs,outputDirPath,"2019",COL_NUMS_100KB,epsilon);

    //5min
//    processMultipleTracesAndGetAvgMetrics<ROW_NUM,BASE_BIT_NUM,6>(traceFileDir2019,fileNum,step,pktsOfFiles,actualFlowSpreadsOfFiles,actualItemNumOfFiles,topKs,outputDirPath,"2019",COL_NUMS_300KB,epsilon);
//    processMultipleTracesAndGetAvgMetrics<ROW_NUM,BASE_BIT_NUM,6>(traceFileDir2019,fileNum,step,pktsOfFiles,actualFlowSpreadsOfFiles,actualItemNumOfFiles,topKs,outputDirPath,"2019",COL_NUMS_250KB,epsilon);
//    processMultipleTracesAndGetAvgMetrics<ROW_NUM,BASE_BIT_NUM,7>(traceFileDir2019,fileNum,step,pktsOfFiles,actualFlowSpreadsOfFiles,actualItemNumOfFiles,topKs,outputDirPath,"2019",COL_NUMS_200KB,epsilon);
//    processMultipleTracesAndGetAvgMetrics<ROW_NUM,BASE_BIT_NUM,8>(traceFileDir2019,fileNum,step,pktsOfFiles,actualFlowSpreadsOfFiles,actualItemNumOfFiles,topKs,outputDirPath,"2019",COL_NUMS_150KB,epsilon);
//    processMultipleTracesAndGetAvgMetrics<ROW_NUM,BASE_BIT_NUM,9>(traceFileDir2019,fileNum,step,pktsOfFiles,actualFlowSpreadsOfFiles,actualItemNumOfFiles,topKs,outputDirPath,"2019",COL_NUMS_100KB,epsilon);
    return 0;
}
