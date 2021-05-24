#include "LDTree.h" 

LDTree::LDTree(const Snp& snp, const Location& location)
{
    nodes_.push_back(snp);
    genomeLocations_.push_back(location);
}

LDTree::LDTree(const LDTree& cpy)
{
    nodes_ = cpy.nodes_;
    genomeLocations_ = cpy.genomeLocations_;
}

LDTree::LDTree(const LDTree & t1, const LDTree & t2)
{
    nodes_.reserve(1 + t1.size() + t2.size());
    genomeLocations_.reserve(1 + t1.size() + t2.size());
    
    //Create new root node
    nodes_.push_back(Snp(t1.nodes_[0], t2.nodes_[0]));
    genomeLocations_.push_back(Location(t1.genomeLocations_[0], t2.genomeLocations_[0]));
    
    //Push back the sub trees
    
    for(int i=0; i<floor(log2(t1.size())+1); ++i)
    {
        for(int j=(1<<(i))-1; j< (1<<(i+1))-1; ++j)
        {
           // cout<<j<<endl;
            nodes_.push_back(t1.nodes_[j]);
            genomeLocations_.push_back(t1.genomeLocations_[j]);
        }
        for(int j=(1<<(i))-1; j< (1<<(i+1))-1; ++j)
        {
            nodes_.push_back(t2.nodes_[j]);
            genomeLocations_.push_back(t2.genomeLocations_[j]);
        }
    }
}

bool LDTree::empty()const
{
    return nodes_.empty();
}

size_t LDTree::size()const
{
    return nodes_.size();
}

const Snp & LDTree::getRoot()const
{
    return nodes_[0];
}

int LDTree::computeDifferences(const LDTree & other)const
{
    return nodes_[0].computeDifferences(other.nodes_[0]);
}

void LDTree::clear()
{
    nodes_.clear();
}

void LDTree::epistasisTest(const LDTree & other, TopSnpList & topSnpList)const
{
    //Vector to use as a stack for tests
    std::vector<std::pair<int, int> > s;
    
    size_t leafStart1 = nodes_.size()/2;
    size_t leafStart2 = other.nodes_.size()/2;

    uint64_t  localInternalTestsDone = 0;
    uint64_t  localLeaftTestsDone = 0;
    
    s.push_back(std::make_pair(0,0));
    while(!s.empty())
    {
        float cutOff=topSnpList.getCutoff();

        std::pair<int,int> c = s.back();
        s.pop_back();

        //cout<<c.first<<" "<<c.second<<endl;
        
        float score = nodes_[c.first].epistasisTest(other.nodes_[c.second]);

        int l1= (c.first<<1) +1;
        int r1= (c.first<<1) +2;
        int l2= (c.second<<1) +1;
        int r2= (c.second<<1) +2;

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
           // ++other.internalTestsDone_;
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
                topSnpList.attemptInsert(nodes_[c.first].getIndex(), other.nodes_[c.second].getIndex(), score);
                localLeaftTestsDone += 2;
            }
                
        }
        //If score below cutoff
            //Do nothing
    }
    topSnpList.incrementInternalTestsCounter(localInternalTestsDone);
    topSnpList.incrementLeafTestsCounter(localLeaftTestsDone);
}

void LDTree::epistasisTestNoTrees(const LDTree & other, TopSnpList & topSnpList)const
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
