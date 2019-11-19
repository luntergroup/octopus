// Copyright (c) 2015-2019 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef phylogeny_hpp
#define phylogeny_hpp

#include <unordered_map>
#include <memory>
#include <cstddef>
#include <stack>
#include <type_traits>

#include <boost/optional.hpp>

namespace octopus {

template <typename Label, typename T = boost::optional<int>>
class Phylogeny
{
public:
    using LabelType = Label;
    using ValueType = T;
    
    struct Group
    {
        Label id;
        T value = boost::none;
    };
    
    Phylogeny() = default;
    
    Phylogeny(Group founder);
    
    Phylogeny(const Phylogeny&);
    Phylogeny& operator=(Phylogeny);
    Phylogeny(Phylogeny&&) = default;
    
    ~Phylogeny() = default;
    
    std::size_t size() const noexcept;
    bool empty() const noexcept;
    
    void clear() noexcept;
    void clear(const Label& id);
    
    std::vector<Group> groups() const;
    
    Group& set_founder(Group founder);
    Group& add_descendant(Group group, const Label& ancestor_id);
    
    const Group& group(const Label& id) const;
    Group& group(const Label& id);
    
    const Group& founder() const noexcept;
    const Group& ancestor(const Label& id) const;
    
    unsigned num_descendants(const Label& id) const noexcept;
    const Group& descendant1(const Label& id) const;
    const Group& descendant2(const Label& id) const;
    
    template <typename UnaryFunction>
    Phylogeny<Label, std::result_of_t<UnaryFunction(T)>>
    transform(UnaryFunction op) const;

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
Phylogeny<Label, T>::Phylogeny(Group founder)
{
    set_founder(std::move(founder));
}

template <typename Label, typename T>
Phylogeny<Label, T>::Phylogeny(const Phylogeny& other)
{
    if (other.tree_) {
        std::stack<TreeNode*> to_visit {};
        to_visit.push(other.tree_->descendant2.get());
        to_visit.push(other.tree_->descendant1.get());
        set_founder(other.tree_->group);
        while (!to_visit.empty()) {
            if (to_visit.top() != nullptr) {
                const TreeNode* visted {to_visit.top()};
                to_visit.pop();
                add_descendant(visted->group, visted->ancestor->group.id);
                to_visit.push(visted->descendant2.get());
                to_visit.push(visted->descendant1.get());
            } else {
                to_visit.pop();
            }
        }
    }
}

template <typename Label, typename T>
Phylogeny<Label, T>& Phylogeny<Label, T>::operator=(Phylogeny other)
{
    std::swap(*this, other);
    return *this;
}

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
std::vector<typename Phylogeny<Label, T>::Group> Phylogeny<Label, T>::groups() const
{
    std::vector<Group> result {};
    result.reserve(size());
    for (const auto& p : nodes_) {
        result.push_back(p.second->group);
    }
    return result;
}

template <typename Label, typename T>
typename Phylogeny<Label, T>::Group& Phylogeny<Label, T>::set_founder(Group founder)
{
    if (!tree_) {
        TreeNode node {std::move(founder), nullptr, nullptr, nullptr};
        tree_ = std::make_unique<TreeNode>(std::move(node));
    } else {
        nodes_.erase(founder.id);
        tree_->group = std::move(founder);
    }
    return nodes_.emplace(founder.id, std::addressof(*tree_)).first->second->group;
}

template <typename Label, typename T>
typename Phylogeny<Label, T>::Group& Phylogeny<Label, T>::add_descendant(Group group, const Label& ancestor_id)
{
    TreeNode* ancestor {nodes_.at(ancestor_id)};
    if (ancestor->descendant1) {
        if (ancestor->descendant2) {
            throw std::runtime_error {"Cannot add descendant to node with two descendants"};
        }
        TreeNode node {std::move(group), ancestor, nullptr, nullptr};
        ancestor->descendant2 = std::make_unique<TreeNode>(std::move(node));
        return nodes_.emplace(ancestor->descendant2->group.id, std::addressof(*ancestor->descendant2)).first->second->group;
    } else {
        TreeNode node {std::move(group), ancestor, nullptr, nullptr};
        ancestor->descendant1 = std::make_unique<TreeNode>(std::move(node));
        return nodes_.emplace(ancestor->descendant1->group.id, std::addressof(*ancestor->descendant1)).first->second->group;
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
    return tree_->group;
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

template <typename Label, typename T>
template <typename UnaryFunction>
Phylogeny<Label, std::result_of_t<UnaryFunction(T)>>
Phylogeny<Label, T>::transform(UnaryFunction op) const
{
    using T2 = std::result_of_t<UnaryFunction(T)>;
    Phylogeny<Label, T2> result {};
    std::stack<TreeNode*> to_visit {};
    to_visit.push(this->tree_->descendant2.get());
    to_visit.push(this->tree_->descendant1.get());
    const auto transformer = [&] (const Group& old) -> typename Phylogeny<Label, T2>::Group { return {old.id, op(old.value)}; };
    result.set_founder(transformer(this->tree_->group));
    while (!to_visit.empty()) {
        if (to_visit.top() != nullptr) {
            const TreeNode* visted {to_visit.top()};
            to_visit.pop();
            result.add_descendant(transformer(visted->group), visted->ancestor->group.id);
            to_visit.push(visted->descendant2.get());
            to_visit.push(visted->descendant1.get());
        } else {
            to_visit.pop();
        }
    }
    return result;
}

} // namespace octopus

#endif
