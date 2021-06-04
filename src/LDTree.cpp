#include "LDTree.h" 

LDTree::LDTree(const Snp& snp, const Location& location){
    nodes_.push_back(snp);
    genomeLocations_.push_back(location);
    topSnpList_ = nullptr;//topSnpList;
}

LDTree::LDTree(const LDTree& cpy){
    nodes_ = cpy.nodes_;
    genomeLocations_ = cpy.genomeLocations_;
    topSnpList_ = cpy.topSnpList_;
}

LDTree::LDTree(const LDTree & t1, const LDTree & t2){
    nodes_.reserve(1 + t1.size() + t2.size());
    genomeLocations_.reserve(1 + t1.size() + t2.size());
    topSnpList_ = t1.topSnpList_;
    
    //Create new root node
    nodes_.push_back(Snp(t1.nodes_[0], t2.nodes_[0]));
    genomeLocations_.push_back(Location());

    //Push back the sub trees
    for(int i=0; i<floor(log2(t1.size())+1); ++i){
        for(int j=(1<<(i))-1; j< (1<<(i+1))-1; ++j){
            nodes_.push_back(t1.nodes_[j]);
            genomeLocations_.push_back(t1.genomeLocations_[j]);
        }
        for(int j=(1<<(i))-1; j< (1<<(i+1))-1; ++j){
            nodes_.push_back(t2.nodes_[j]);
            genomeLocations_.push_back(t2.genomeLocations_[j]);
        }
    }
}

bool LDTree::empty()const{
    return nodes_.empty();
}

size_t LDTree::size()const{
    return nodes_.size();
}

const Snp & LDTree::getRoot()const{
    return nodes_[0];
}

ID_Sample LDTree::computeDifferences(const LDTree & other)const{
    return nodes_[0].computeDifferences(other.nodes_[0]);
}

bool LDTree::validMerge(const LDTree& other, ID_Sample maxDiff)const{
    //if ((!ldtrees_[j].empty()) && (ldtrees_[i].size() == ldtrees_[j].size())) { //check on equivalent tree sizes
      //  if (ldtrees_[i].computeDifferences(ldtrees_[j]) < allowedDifferences) { //check on tree overlap
    return (!empty() && size() == other.size()) && computeDifferences(other) < maxDiff;
}

void LDTree::clear(){
    nodes_.clear();
}

void LDTree::epistasisTest(const LDTree & other)const{
    //Vector to use as a stack for tests
    std::vector<std::pair<size_t, size_t> > s;
    
    size_t leafStart1 = nodes_.size()/2;
    size_t leafStart2 = other.nodes_.size()/2;

    uint64_t  localInternalTestsDone = 0;
    uint64_t  localLeaftTestsDone = 0;
    
    CTable2 cTable;

    s.push_back(std::make_pair(0,0));
    while(!s.empty())
    {
        float cutOff=topSnpList_->getCutoff();

        std::pair<size_t, size_t> c = s.back();
        s.pop_back();
      
       // float score = nodes_[c.first].epistasisTest(other.nodes_[c.second]);
        Snp::fillTable(cTable, nodes_[c.first], other.nodes_[c.second]);
        float score = cTable.chi2();

        size_t l1= (c.first<<1) +1;
        size_t r1= (c.first<<1) +2;
        size_t l2= (c.second<<1) +1;
        size_t r2= (c.second<<1) +2;

        //If both not at both leaves
        if(c.first<leafStart1 && c.second<leafStart2)
        {
            if(score >= cutOff)
            {
                s.push_back(std::make_pair(l1, l2));
                s.push_back(std::make_pair(l1, r2));
                s.push_back(std::make_pair(r1, l2));
                s.push_back(std::make_pair(r1, r2));
            }
            localInternalTestsDone += 2;
        }
        else if(c.first<leafStart1)
        {
            if(score >= cutOff)
            {
                s.push_back(std::make_pair(l1, c.second ));
                s.push_back(std::make_pair(r1, c.second));
            }
            ++localInternalTestsDone;
            ++localLeaftTestsDone;
        }
        else if(c.second<leafStart2)
        {
            if(score >= cutOff)
            {
                s.push_back(std::make_pair(c.first, l2));
                s.push_back(std::make_pair(c.first, r2));
            }    
            ++localLeaftTestsDone;
            ++localInternalTestsDone;
        }
        //If at both leaves
        else
        {
            //check to make sure not estimated as being in LD
            if(!genomeLocations_[c.first].inLinkageDisequilibrium(other.genomeLocations_[c.second]) )
            {
                topSnpList_->attemptInsert(nodes_[c.first].getIndex(), other.nodes_[c.second].getIndex(), score);
                localLeaftTestsDone += 2;
            }
                
        }
        //If score below cutoff
            //Do nothing
    }
    topSnpList_->incrementTestCounter(TestCounter{ localInternalTestsDone,  localLeaftTestsDone });
}

/*void LDTree::epistasisTestNoTrees(const LDTree& other, TopSnpList& topSnpList)const
{
   // cout<<"HERE"<<endl;
    size_t leafStart1 = nodes_.size()/2;
    size_t leafStart2 = other.nodes_.size()/2;
    
    int n1 = rand() % (leafStart1+1) + leafStart1;
    int n2 = rand() % (leafStart2+1) + leafStart2;
   // cout<<n1<<" "<<n2<<endl;
    if(!genomeLocations_[n1].inLinkageDisequilibrium(other.genomeLocations_[n2]) )
    {
        float score = nodes_[n1].epistasisTest(other.nodes_[n2]);
        topSnpList.attemptInsert(nodes_[n1].getIndex(), other.nodes_[n2].getIndex(), score);
    }
    topSnpList.incrementLeafTestsCounter(2);
}
*/