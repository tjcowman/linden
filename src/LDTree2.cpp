#include "LDTree2.h" 

LDTree::LDTree(const Snp& snp, const Location& location) : snps_(Graph<Snp, ID_Snp>(snp)), locations_(std::vector<Location>{location}) {
    topSnpList_ = nullptr;
    root_ = 0;
}

LDTree::LDTree( LDTree& t1,  LDTree& t2)  {
    root_ = 0;
   
    //Create new representative snp
    Snp newRoot(t1.getRoot(), t2.getRoot());

    //Update the locations data
    locations_ = std::vector<Location>();
    locations_.insert(locations_.end(), std::make_move_iterator(t1.locations_.begin()), std::make_move_iterator(t1.locations_.end()));
    locations_.insert(locations_.end(), std::make_move_iterator(t2.locations_.begin()), std::make_move_iterator(t2.locations_.end()));
   
    //Update the Snp Tree (Graph) structure
    snps_ = Graph<Snp, ID_Snp>::joinToRoot(newRoot, t1.snps_, t2.snps_);

    //set the correct topSnpList pointer fail with error if the subtrees point to differnt top lists
    if (t1.topSnpList_ != t2.topSnpList_) {
        std::cerr << "subtrees point to different top snp lists during merge" << std::endl;
        exit(1);
    }
    else {
        topSnpList_ = t1.topSnpList_;
    }


}

bool LDTree::empty()const{
    return snps_.empty();
}

size_t LDTree::size()const {
    return snps_.size();
}

const Snp& LDTree::getRoot()const {
    return snps_.getElement(0);
}

std::vector<ID_Snp> LDTree::getChildren(ID_Snp i)const {
    return snps_.getOutgoingIds(i);
}

bool LDTree::isLeaf(ID_Snp i)const {
    return snps_.isTerminal(i);
}

const Snp& LDTree::getSnp(ID_Snp i)const {
    return snps_.getElement(i);
}


ID_Sample LDTree::computeDifferences(const LDTree& other)const {
    return getRoot().computeDifferences(other.getRoot());
}

bool LDTree::validMerge(const LDTree& other, ID_Sample maxDiff)const {
    //The restriction used to be to enusre binary trees, with this DS it is no longer required. However a useful hueristic should be investigated.
    if (size() == 0 || other.size() == 0) return false;
    else
        return computeDifferences(other) < maxDiff;
}


void LDTree::clear() {
    root_ = ID_Invalid::Snp;
    snps_.clear();
    topSnpList_ = nullptr;
}


void LDTree::epistasisTest(const LDTree& other)const {
    uint64_t localInternalTestsDone = 0;
    uint64_t localLeaftTestsDone = 0;
    CTable2 cTable;

    //Vector to use as a stack for tests
    std::vector<std::pair<ID_Snp, ID_Snp> > s;
    s.push_back({ 0, 0 });

    while (!s.empty()) {
        float cutOff = topSnpList_->getCutoff();      
        auto c = s.back();

        s.pop_back();


        Snp::fillTable(cTable, snps_.getElement(c.first), other.snps_.getElement(c.second));
        float score = cTable.chi2();
     
        auto lChildren = getChildren(c.first);
        auto rChildren = other.getChildren(c.second);
        

        //If both not at leaves
        if(!isLeaf(c.first) && !other.isLeaf(c.second)) {
            if (score >= cutOff) {
                s.push_back(std::make_pair(lChildren[0], rChildren[0]));
                s.push_back(std::make_pair(lChildren[0], rChildren[1]));
                s.push_back(std::make_pair(lChildren[1], rChildren[0]));
                s.push_back(std::make_pair(lChildren[1], rChildren[1]));
            }
            localInternalTestsDone += 2;
        }
        else if(!isLeaf(c.first)){
            if (score >= cutOff){
                s.push_back(std::make_pair(lChildren[0], c.second));
                s.push_back(std::make_pair(lChildren[1], c.second));

            }
            ++localInternalTestsDone;
            ++localLeaftTestsDone;
        }
        else if(!other.isLeaf(c.second)){
            if (score >= cutOff){
                s.push_back(std::make_pair(c.first, rChildren[0]));
                s.push_back(std::make_pair(c.first, rChildren[1]));
            }
            ++localLeaftTestsDone;
            ++localInternalTestsDone;
        }
        //If at both leaves
        else{
            //check to make sure not estimated as being in LD
            if (!locations_[c.first].inLinkageDisequilibrium(other.locations_[c.second])){

                topSnpList_->attemptInsert(snps_.getElement(c.first).getIndex(), other.snps_.getElement(c.second).getIndex(), score);

                localLeaftTestsDone += 2;
           }

        }
        //If score below cutoff
            //Do nothing
    }
    
    topSnpList_->incrementTestCounter(TestCounter{ localInternalTestsDone,  localLeaftTestsDone });
}