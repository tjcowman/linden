#include "LDTree2.h" 

LDTree::LDTree() {

}

LDTree::LDTree(const Snp& snp, const Location& location) : 
    snps_(Graph<Snp, ID_Snp>(snp)), 
    locations_(Graph<Location, ID_Snp>(location)) 
{
    topSnpList_ = nullptr;
    root_ = 0;
}

LDTree::LDTree( LDTree& t1,  LDTree& t2)  {
    root_ = 0;
   
    //Create new representative snp
    Snp newRoot(t1.getRoot(), t2.getRoot());

    //Update the locations data
   // locations_ = std::vector<Location>();
    //locations_.insert(locations_.end(), std::make_move_iterator(t1.locations_.begin()), std::make_move_iterator(t1.locations_.end()));
    //locations_.insert(locations_.end(), std::make_move_iterator(t2.locations_.begin()), std::make_move_iterator(t2.locations_.end()));
   
    locations_ = Graph<Location, ID_Snp>::joinToRoot(Location(), t1.locations_, t2.locations_);

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

bool LDTree::operator==(const LDTree& lhs)const {
    return std::tie(root_, snps_, locations_) == std::tie(lhs.root_, lhs.snps_, lhs.locations_);
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
            //if (!locations_[c.first].inLinkageDisequilibrium(other.locations_[c.second])){
            if (!locations_.getElement(c.first).inLinkageDisequilibrium(other.locations_.getElement(c.second))){

                topSnpList_->attemptInsert(snps_.getElement(c.first).getIndex(), other.snps_.getElement(c.second).getIndex(), score);

                localLeaftTestsDone += 2;
           }

        }
        //If score below cutoff
            //Do nothing
    }
    
    topSnpList_->incrementTestCounter(TestCounter{ localInternalTestsDone,  localLeaftTestsDone });
}

void LDTree::to_serial(std::ostream& os, const LDTree& e){
 /*  
    //Root node
    os.write(reinterpret_cast<const char*>(&e.root_), sizeof(ID_Snp));
   // Graph<Snp, ID_Snp> snps_;
    Graph<Snp, ID_Snp>::to_serial(os, e.snps_);

   // std::vector<Location> locations_;
    ID_Snp num_locations = e.locations_.size();
    os.write(reinterpret_cast<const char*>(&num_locations), sizeof(ID_Snp));
    os.write(reinterpret_cast<const char*>(&e.locations_[0]), num_locations*sizeof(Location));
*/

}

LDTree LDTree::from_serial(std::istream& is) {
 /*   LDTree e;

    is.read(reinterpret_cast<char*>(&e.root_), sizeof(ID_Snp));
    e.snps_ = Graph<Snp, ID_Snp>::from_serial(is);

    ID_Snp num_locations;
    is.read(reinterpret_cast<char*>(&num_locations), sizeof(ID_Snp));

    e.locations_.resize(num_locations);
    is.read(reinterpret_cast<char*>(&e.locations_[0]), num_locations*sizeof(Location));

    return e;*/
}