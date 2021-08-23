#include "LDTree2.h" 

LDTree::LDTree(const Snp& snp, const Location& location) : snps_(Graph<Snp, ID_Snp>(snp)), locations_(std::vector<Location>{location}) {
    topSnpList_ = nullptr;
    root_ = 0;
}

//TODO initlist tmp to compile
LDTree::LDTree( LDTree& t1,  LDTree& t2)  {
    //TODO
    root_ = 0;

    //Create new representative snp
    Snp newRoot(t1.getRoot(), t2.getRoot());

    //Update the locations data

    //Update the Snp Tree (Graph) structure
    snps_ = Graph<Snp, ID_Snp>::joinToRoot(newRoot, t1.snps_, t2.snps_);


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

const Snp& LDTree::getLeft(ID_Snp id)const {

}

const Snp& LDTree::getRight(ID_Snp id)const {

}

ID_Sample LDTree::computeDifferences(const LDTree& other)const {
    return getRoot().computeDifferences(other.getRoot());
}

bool LDTree::validMerge(const LDTree& other, ID_Sample maxDiff)const {
    //TODO
    //The restriction used to be to enusre binary trees, with this DS it is no longer required. However a useful hueristic should be investigated.
    if (size() == other.size())
        return true;
    else
        return false;
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
/*
    //Vector to use as a stack for tests
    std::vector<std::pair<size_t, size_t> > s;

    size_t leafStart1 = nodes_.size() / 2;
    size_t leafStart2 = other.nodes_.size() / 2;



    s.push_back(std::make_pair(0, 0));
    while (!s.empty())
    {
        float cutOff = topSnpList_->getCutoff();

        std::pair<size_t, size_t> c = s.back();
        s.pop_back();

        // float score = nodes_[c.first].epistasisTest(other.nodes_[c.second]);
        Snp::fillTable(cTable, nodes_[c.first], other.nodes_[c.second]);
        float score = cTable.chi2();

        size_t l1 = (c.first << 1) + 1;
        size_t r1 = (c.first << 1) + 2;
        size_t l2 = (c.second << 1) + 1;
        size_t r2 = (c.second << 1) + 2;

        //If both not at both leaves
        if (c.first < leafStart1 && c.second < leafStart2)
        {
            if (score >= cutOff)
            {
                s.push_back(std::make_pair(l1, l2));
                s.push_back(std::make_pair(l1, r2));
                s.push_back(std::make_pair(r1, l2));
                s.push_back(std::make_pair(r1, r2));
            }
            localInternalTestsDone += 2;
        }
        else if (c.first < leafStart1)
        {
            if (score >= cutOff)
            {
                s.push_back(std::make_pair(l1, c.second));
                s.push_back(std::make_pair(r1, c.second));
            }
            ++localInternalTestsDone;
            ++localLeaftTestsDone;
        }
        else if (c.second < leafStart2)
        {
            if (score >= cutOff)
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
            if (!genomeLocations_[c.first].inLinkageDisequilibrium(other.genomeLocations_[c.second]))
            {
                topSnpList_->attemptInsert(nodes_[c.first].getIndex(), other.nodes_[c.second].getIndex(), score);
                localLeaftTestsDone += 2;
            }

        }
        //If score below cutoff
            //Do nothing
    }
    */
    topSnpList_->incrementTestCounter(TestCounter{ localInternalTestsDone,  localLeaftTestsDone });
}