
#include "LDTree.h"

namespace Linden::Core
{
    const Snp& LDTree::getRoot() const {
        return snps_.getElement(0);
    }

    bool LDTree::isLeaf(Genetics::Id::Snp i) const {
        return snps_.isTerminal(i);
    }

    Genetics::Id::Sample LDTree::computeDifferences(const LDTree& other) const {
        return getRoot().computeDifferences(other.getRoot());
    }

    bool LDTree::validMerge(const LDTree& other, Genetics::Id::Sample maxDiff) const {
        //The restriction used to be to enusre binary trees, with this DS it is no longer required. However a useful hueristic should be investigated.
        if (size() == 0 || other.size() == 0) return false;
        else
            return computeDifferences(other) < maxDiff;
    }

    void LDTree::epistasisTest(const LDTree& other) const {
        uint64_t localInternalTestsDone = 0;
        uint64_t localLeaftTestsDone = 0;
        Statistics::ContingencyTable<9> cTable;

        //Vector to use as a stack for tests
        std::vector<std::pair<Genetics::Id::Snp, Genetics::Id::Snp> > s;
        s.push_back({ 0, 0 });

        while (!s.empty()) {
            float cutOff = topSnpList_->getCutoff();
            auto c = s.back();

            s.pop_back();


            Snp::fillTable(cTable, snps_.getElement(c.first), other.snps_.getElement(c.second));
            float score = cTable.Chi2();

            //auto lChildren = getChildren(c.first);
            //auto rChildren = other.getChildren(c.second);

            auto lChildren = getChildren2(c.first);
            auto rChildren = other.getChildren2(c.second);

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
    // return LDTree();
    }
}