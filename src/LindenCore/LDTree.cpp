#include "LDTree.h"



namespace Linden::Core
{
        
    LDTree::LDTree() {

    }

    /*
    LDTree::LDTree(const LDTree& cpy){
        nodes_ = cpy.nodes_;
        genomeLocations_ = cpy.genomeLocations_;
        topSnpList_ = cpy.topSnpList_;
    }
    */

    bool LDTree::operator==(const LDTree& lhs) const {
        return std::tie(nodes_, genomeLocations_) == std::tie(lhs.nodes_, lhs.genomeLocations_);
    }

    LDTree::LDTree( LDTree & t1,  LDTree & t2){
        nodes_.reserve(1 + t1.size() + t2.size());
        genomeLocations_.reserve(1 + t1.size() + t2.size());
        topSnpList_ = t1.topSnpList_;

        //Create new root node
        nodes_.push_back(Snp(t1.nodes_[0], t2.nodes_[0]));
        genomeLocations_.push_back(Genetics::Location());

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

    bool LDTree::empty() const{
        return nodes_.empty();
    }

    size_t LDTree::size() const{
        return nodes_.size();
    }

    const Snp & LDTree::getRoot() const{
        return nodes_[0];
    }

    Genetics::Id::Sample LDTree::computeDifferences(const LDTree & other) const{
        return nodes_[0].computeDifferences(other.nodes_[0]);
    }

    //check on equivalent tree sizes and check on tree overlap
    bool LDTree::validMerge(const LDTree& other, Genetics::Id::Sample maxDiff) const{
        return (!empty() && size() == other.size()) && computeDifferences(other) < maxDiff;
    }

    void LDTree::clear(){
        nodes_.clear();
    }

    void LDTree::epistasisTest(const LDTree & other) const{
        //Vector to use as a stack for tests
        std::vector<std::pair<size_t, size_t> > s;

        size_t leafStart1 = nodes_.size()/2;
        size_t leafStart2 = other.nodes_.size()/2;

        uint64_t  localInternalTestsDone = 0;
        uint64_t  localLeaftTestsDone = 0;

        ContingencyTable<9> cTable;

        s.push_back(std::make_pair(0,0));
        while(!s.empty())
        {
            float cutOff=topSnpList_->getCutoff();

            std::pair<size_t, size_t> c = s.back();
            s.pop_back();

        // float score = nodes_[c.first].epistasisTest(other.nodes_[c.second]);
            Snp::fillTable(cTable, nodes_[c.first], other.nodes_[c.second]);
            float score = cTable.Chi2();

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

    /*void LDTree::epistasisTestNoTrees(const LDTree& other, TopSnpList& topSnpList) const
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

    void LDTree::to_serial(std::ostream& os, const LDTree& e){
    /*
        //Root node
        os.write(reinterpret_cast<const char*>(&e.root_), sizeof(Genetics::Id::Snp));
    // Graph<Snp, Genetics::Id::Snp> snps_;
        Graph<Snp, Genetics::Id::Snp>::to_serial(os, e.snps_);

    // std::vector<Location> locations_;
        Genetics::Id::Snp num_locations = e.locations_.size();
        os.write(reinterpret_cast<const char*>(&num_locations), sizeof(Genetics::Id::Snp));
        os.write(reinterpret_cast<const char*>(&e.locations_[0]), num_locations*sizeof(Location));
    */

    }

    LDTree LDTree::from_serial(std::istream& is) {
        return LDTree(Snp(0, std::vector<Bitwise::Genotype>()),Genetics::Location());
    /*   LDTree e;

        is.read(reinterpret_cast<char*>(&e.root_), sizeof(Genetics::Id::Snp));
        e.snps_ = Graph<Snp, Genetics::Id::Snp>::from_serial(is);

        Genetics::Id::Snp num_locations;
        is.read(reinterpret_cast<char*>(&num_locations), sizeof(Genetics::Id::Snp));

        e.locations_.resize(num_locations);
        is.read(reinterpret_cast<char*>(&e.locations_[0]), num_locations*sizeof(Location));

        return e;*/
    }
}
