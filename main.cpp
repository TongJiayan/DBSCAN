//
// Created by JoeyTong on 2021/5/17.
//

#define ATTR_SIZE 4

#define EPS 0.135 //0.135 0.183,0.228
#define MINPOINTS 5

#include <iostream>
#include <fstream>
#include <random>
#include <sstream>
#include <map>
#include <vector>
#include <cmath>
using namespace std;

struct Data{
    vector<double> attr;
    int label;
};

/* 读取数据文件存到totalData */
void readData(const string& filename, vector<Data>& totalData){
    ifstream fin(filename, ios::in);
    if(fin.fail()){
        cout<<"打开数据文件失败！"<<endl;
        exit(1);
    }
    string line;
    while(getline(fin, line)){
        istringstream sin(line);
        string field;
        Data data;
        data.attr.clear();
        while (getline(sin, field, ',')){
            if(data.attr.size() < ATTR_SIZE)
                data.attr.push_back(stod(field));
            else
                data.label = stoi(field);
        }
        totalData.push_back(data);
    }
    shuffle(totalData.begin(), totalData.end(), std::mt19937(std::random_device()()));
}

/* 计算两个样本之间的欧氏距离 */
double computeDistance(Data data1, Data data2){
    double tmp=0;
    for(int i=0;i<ATTR_SIZE;i++)
        tmp += pow((data1.attr[i]-data2.attr[i]),2);
    return sqrt(tmp);
}

void getNeighborSet(vector<Data> totalData, vector<vector<int>>& neighborSet){
    vector<int> neighbor;
    double distance = 0;
    for(int i=0;i<totalData.size();i++){
        for(int j=0;j<totalData.size();j++){
            distance= computeDistance(totalData[i],
                                      totalData[j]);
            if(distance<=EPS && distance!=0)           //<= or <
                neighbor.push_back(j);
        }
        neighborSet.push_back(neighbor);
        neighbor.clear();
    }
}

/* 邻域样本数超过MINPOITS的样本归入核心对象集 */
void getKernelSet(vector<vector<int>> neighborSet, vector<int>& kernelSet){
    for(int i=0;i<neighborSet.size();i++)
        if(neighborSet[i].size()>=MINPOINTS)
            kernelSet.push_back(i);
}

void cluster(vector<vector<int>> neighborSet, vector<int> kernelSet, vector<vector<int>>& clusters){
    vector<int> unClusteredData;
    for(int i=0;i<neighborSet.size();i++)
        unClusteredData.push_back(i);
    int currentKernel=0;
    vector<int> clusterQueue;
    vector<int> cluster;
    while(!kernelSet.empty()){
        // 备份当前未被访问的样本集合（在totalData中位置值）
        vector<int> unClusteredDataOld;
        unClusteredDataOld.assign(unClusteredData.begin(),unClusteredData.end());
        // 选取一个核心对象（在totalData中位置值）塞入Q中
        currentKernel=kernelSet.front();
        clusterQueue.push_back(currentKernel);
        // 从未被访问的样本集合中删除该核心对象对应的样本
        for(int i=0;i<unClusteredData.size();i++)
            if(unClusteredData[i]==currentKernel){
                unClusteredData.erase(unClusteredData.begin()+i);
                break;
            }
        while(!clusterQueue.empty()){
            int sampleIndex = clusterQueue.back();
            clusterQueue.pop_back();
            if(neighborSet[sampleIndex].size()>=MINPOINTS){ //如果是核心对象，将其邻域内的样本归入该聚类蔟
                for(auto neighborPoint : neighborSet[sampleIndex]){
                    auto it = find(unClusteredData.begin(),unClusteredData.end(),neighborPoint);
                    if(it != unClusteredData.end()){
                        clusterQueue.push_back(neighborPoint);
                        unClusteredData.erase(it);
                    }
                }
            }
        }
        // 将unClusteredDataOld和unClusteredData集合的差集推入cluster，称为一个簇
        for(int i=0;i<unClusteredDataOld.size();i++){
            if(find(unClusteredData.begin(),unClusteredData.end(),unClusteredDataOld[i])==unClusteredData.end())
                cluster.push_back(unClusteredDataOld[i]);
        }
        // 更新核心集，将已聚类成簇的样本移出核心集
        for(int i=kernelSet.size()-1;i>=0;i--){
            for(int j=0;j<cluster.size();j++){
                if(cluster[j]==kernelSet[i])
                {
                    kernelSet.erase(kernelSet.begin()+i);
                    break;
                }

            }
        }
        clusters.push_back(cluster);
        cluster.clear();
        unClusteredDataOld.clear();
        clusterQueue.clear();
    }
}

/* DBI */
double evaluation(vector<vector<int>> clusters, vector<Data> totalData)
{
    vector<double> avgInterDist;    //簇内样本间的平均距离
    vector<double> centerPoint(ATTR_SIZE,0);
    vector<vector<double>> centerPoints(clusters.size(),centerPoint); // 簇中心点
    // 计算簇内样本间的平均距离 和 各簇的中心点
    for(int i=0;i<clusters.size();i++)
    {
        double dist=0;
        int clusterSize=clusters[i].size();
        for(auto dataIndex1 : clusters[i])
        {
            for(auto dataIndex2 : clusters[i])
                dist += computeDistance(totalData[dataIndex1], totalData[dataIndex2]);
            for(int j=0;j<ATTR_SIZE;j++)
                centerPoints[i][j] += totalData[dataIndex1].attr[j];
        }
        avgInterDist.push_back(2*dist/(clusterSize*(clusterSize-1)));
        for(int j=0;j<ATTR_SIZE;j++)
            centerPoints[i][j] /= clusterSize;
    }

    // 计算簇中心点之间的距离，因矩阵对称性，故仅需计算一半，即为三角矩阵
    vector<double> centerDist(clusters.size(),0);
    vector<vector<double>> centerDists(clusters.size(),centerDist);
    for(int i=0;i<centerPoints.size();i++)
    {
        for(int j=i+1;j<centerPoints.size();j++)
        {
            for(int k=0;k<ATTR_SIZE;k++)
                centerDists[i][j] += pow(centerPoints[i][k]-centerPoints[j][k],2);
            centerDists[i][j] = sqrt(centerDists[i][j]);
        }
    }

    // 按照DBI计算公式计算
    double DBI=0;
    for(int i=0;i<clusters.size();i++)
    {
        double max = -1;
        for(int j=0;j<clusters.size();j++)
        {
            if(j==i)
                continue;
            double tmp = (avgInterDist[i]+avgInterDist[j])/(centerDists[i][j]+centerDists[j][i]);
            if(tmp>max)
                max = tmp;
        }
        DBI += max;
    }
    DBI = DBI / clusters.size();
    return DBI;
}

void printData(vector<Data> totalData){
    for(auto data : totalData)
    {
        for(auto attr : data.attr)
            cout<<attr<<" ";
        cout<<"------"<<data.label<<endl;
    }
}

void printCluster(vector<vector<int>> clusters, vector<Data> totalData){
    vector<int> count(clusters.size(),0);
    for(int i=0;i<clusters.size();i++){
        cout<<"CLUSTER-"<<i<<"----------"<<endl;
        for(int j=0;j<clusters[i].size();j++){
            int index=clusters[i][j];
            count[i]++;
            cout<<totalData[index].label<<" ";
            if(count[i]%20==0)
                cout<<endl;
        }
        cout<<endl;
    }
}
int main() {
    const string filename = "iris.csv";
    vector<Data> totalData;
    readData(filename, totalData);

    vector<vector<int>> neighborSet;
    vector<int> kernelSet;
    vector<vector<int>> clusters;
    double DBI=0;

    getNeighborSet(totalData,neighborSet);
    getKernelSet(neighborSet,kernelSet);
    cluster(neighborSet,kernelSet,clusters);

    cout<<"聚类结果："<<endl;
    printCluster(clusters,totalData);
    DBI= evaluation(clusters, totalData);

    cout<<endl<<"DBI="<<DBI<<endl;

    return 0;
}
