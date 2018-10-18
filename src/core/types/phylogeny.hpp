// Copyright (c) 2015-2018 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef phylogeny_hpp
#define phylogeny_hpp

#include <vector>
#include <unordered_map>
#include <memory>
#include <cstddef>

namespace octopus {

template <typename Label, typename T>
class Phylogeny
{
public:
    struct Group
    {
        using MemberArray = std::vector<T>;
        Label id;
        MemberArray members;
    };
    
    Phylogeny() = default;
    
    Phylogeny(Group founder);
    
    Phylogeny(const Phylogeny&)            = default;
    Phylogeny& operator=(const Phylogeny&) = default;
    Phylogeny(Phylogeny&&)                 = default;
    Phylogeny& operator=(Phylogeny&&)      = default;
    
    ~Phylogeny() = default;
    
    std::size_t size() const noexcept;
    bool empty() const noexcept;
    
    void clear() noexcept;
    void clear(const Label& id);
    
    Group& set_founder(Group founder);
    Group& add_descendant(Group group, const Label& ancestor_id);
    
    const Group& group(const Label& id) const;
    Group& group(const Label& id);
    
    const Group& founder() const noexcept;
    const Group& ancestor(const Label& id) const;
    
    unsigned num_descendants(const Label& id) const noexcept;
    const Group& descendant1(const Label& id) const;
    const Group& descendant2(const Label& id) const;

private:
    struct TreeNode
    {
        Group group;
        TreeNode* ancestor = nullptr;
        std::unique_ptr<TreeNode> descendant1 = nullptr, descendant2 = nullptr;
    };
    
    std::unique_ptr<TreeNode> tree_;
    std::unordered_map<Label, TreeNode*> nodes_;
};

template <typename Label, typename T>
std::size_t Phylogeny<Label, T>::size() const noexcept
{
    return nodes_.size();
}

template <typename Label, typename T>
bool Phylogeny<Label, T>::empty() const noexcept
{
    return tree_ == nullptr;
}

template <typename Label, typename T>
void Phylogeny<Label, T>::clear() noexcept
{
    nodes_.clear();
    tree_.reset();
}

template <typename Label, typename T>
void Phylogeny<Label, T>::clear(const Label& id)
{
    if (tree_) {
        if (tree_->group.id == id) {
            clear();
        } else {
            TreeNode* descendant {nodes_.at(id)};
            const TreeNode& ancestor {descendant->ancestor};
            // Remove nodes from cache
            while (descendant != ancestor) {
                if (descendant->descendant1 && nodes_.count(descendant->descendant1->group.id) == 1) {
                    if (descendant->descendant2 && nodes_.count(descendant->descendant2->group.id) == 1) {
                        descendant = descendant->descendant2;
                    } else {
                        descendant = descendant->descendant1;
                    }
                } else {
                    nodes_.erase(descendant->group.id);
                    descendant = descendant->ancestor;
                }
            }
            // Remove descendants from tree
            if (ancestor.descendant1->group.id == id) {
                ancestor.descendant1.reset();
            } else {
                ancestor.descendant2.reset();
            }
        }
    }
}

template <typename Label, typename T>
typename Phylogeny<Label, T>::Group& Phylogeny<Label, T>::set_founder(Group founder)
{
    if (tree_) {
        tree_ = std::make_unique(std::move(founder));
    } else {
        nodes_.erase(founder.id);
        tree_->group = std::move(founder);
    }
    nodes_.emplace(founder.id, std::addressof(*tree_));
}

template <typename Label, typename T>
typename Phylogeny<Label, T>::Group& Phylogeny<Label, T>::add_descendant(Group group, const Label& ancestor_id)
{
    const TreeNode* ancestor {nodes_.at(ancestor_id)};
    if (ancestor->descendant1) {
        ancestor->descendant2 = std::make_unique(std::move(group));
        nodes_.emplace(ancestor->descendant2->group.id, std::addressof(*ancestor->descendant2));
    } else {
        ancestor->descendant1 = std::make_unique(std::move(group));
        nodes_.emplace(ancestor->descendant1->group.id, std::addressof(*ancestor->descendant1));
    }
}

template <typename Label, typename T>
const typename Phylogeny<Label, T>::Group& Phylogeny<Label, T>::group(const Label& id) const
{
    return nodes_.at(id)->group;
}

template <typename Label, typename T>
typename Phylogeny<Label, T>::Group& Phylogeny<Label, T>::group(const Label& id)
{
    return nodes_.at(id)->group;
}

template <typename Label, typename T>
const typename Phylogeny<Label, T>::Group& Phylogeny<Label, T>::founder() const noexcept
{
    tree_->group;
}

template <typename Label, typename T>
const typename Phylogeny<Label, T>::Group& Phylogeny<Label, T>::ancestor(const Label& id) const
{
    return nodes_.at(id)->ancestor->group;
}

template <typename Label, typename T>
unsigned Phylogeny<Label, T>::num_descendants(const Label& id) const noexcept
{
    const TreeNode* ancestor {nodes_.at(id)};
    if (ancestor->descendant2) {
        return 2;
    } else if (ancestor->descendant1) {
        return 1;
    } else {
        return 0;
    }
}

template <typename Label, typename T>
const typename Phylogeny<Label, T>::Group& Phylogeny<Label, T>::descendant1(const Label& id) const
{
    return nodes_.at(id)->descendant1.group;
}

template <typename Label, typename T>
const typename Phylogeny<Label, T>::Group& Phylogeny<Label, T>::descendant2(const Label& id) const
{
    return nodes_.at(id)->descendant2.group;
}

} // namespace octopus

#endif
