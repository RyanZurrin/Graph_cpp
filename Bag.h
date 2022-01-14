//
// Created by Ryan.Zurrin001 on 1/13/2022.
//

#ifndef GRAPH_CPP_BAG_H
#define GRAPH_CPP_BAG_H
#include <bits/stdc++.h>
using namespace std;

template<typename ITEM>
class Node {
public:
    ITEM item;
    Node<ITEM> *next;

    Node(ITEM item_) {
        this->item = item_;
        this->next = nullptr;
    }
};


template <typename ITEM>
class Bag {
    Node<ITEM> *first; // first node in the bag
    int n; // number of items in the bag
public:
    //using ValueType = ITEM;
    Bag() {
        first = nullptr;
        n = 0;
    }
    ~Bag() {
        Node<ITEM> *tmp;
        while (first != nullptr) {
            tmp = first;
            first = first->next;
            delete tmp;
        }
    }
    bool isEmpty() const {
        return first == nullptr;
    }
    int size() const {
        return n;
    }
    void add(const ITEM &item) {
        Node<ITEM> *oldFirst = first;
        first = new Node<ITEM>(item);
        first->item = item;
        first->next = oldFirst;
        n++;
    }
    bool contains(const ITEM &item) const {
        Node<ITEM> *current = first;
        while (current != nullptr) {
            if (current->item == item) {
                return true;
            }
            current = current->next;
        }
        return false;
    }

    void remove(const ITEM &item) {
        Node<ITEM> *current = first;
        Node<ITEM> *previous = nullptr;
        while (current != nullptr) {
            if (current->item == item) {
                if (previous == nullptr) {
                    first = current->next;
                } else {
                    previous->next = current->next;
                }
                delete current;
                n--;
                return;
            }
            previous = current;
            current = current->next;
        }
    }
    void clear() {
        Node<ITEM> *tmp;
        while (first != nullptr) {
            tmp = first;
            first = first->next;
            delete tmp;
        }
        n = 0;
    }
    void print() const {
        Node<ITEM> *current = first;
        while (current != nullptr) {
            std::cout << current->item << " ";
            current = current->next;
        }
        std::cout << std::endl;
    }

    // add an iterator to the bag
    class Iterator {
    public:
        Iterator(Node<ITEM> *current) {
            this->current = current;
        }
        Iterator(const Iterator &other) {
            this->current = other.current;
        }
        Iterator &operator=(const Iterator &other) {
            this->current = other.current;
            return *this;
        }
        Iterator &operator++() {
            this->current = this->current->next;
            return *this;
        }
        Iterator operator++(int) {
            Iterator tmp(*this);
            ++(*this);
            return tmp;
        }
        bool operator==(const Iterator &other) const {
            return this->current == other.current;
        }
        bool operator!=(const Iterator &other) const {
            return this->current != other.current;
        }
        ITEM &operator*() {
            return this->current->item;
        }
        ITEM *operator->() {
            return &(this->current->item);
        }
    private:
        Node<ITEM> *current;
    };

    Iterator begin() {
        return Iterator(first);
    }
    Iterator end() {
        return Iterator(nullptr);
    }
    // overload the [] operator
    ITEM &operator[](int index) {
        Node<ITEM> *current = first;
        for (int i = 0; i < index; i++) {
            current = current->next;
        }
        // make sure the index is valid
        if (current == nullptr) {
            throw std::out_of_range("Index out of range");
        }
        return current->item;
    }
    // overload the () operator
    ITEM &operator()(int index) {
        Node<ITEM> *current = first;
        for (int i = 0; i < index; i++) {
            current = current->next;
        }
        // make sure the index is valid
        if (current == nullptr) {
            throw std::out_of_range("Index out of range");
        }
        return current->item;
    }

    // over load the operator <<
    friend ostream &operator<<(ostream &os, Bag<ITEM> &bag) {
        for (auto it = bag.begin(); it != bag.end(); ++it) {
            os << *it << " ";
        }
        return os;
    }
    // over load the relational operators
    bool operator==(const Bag<ITEM> &other) const {
        if (this->n != other.n) {
            return false;
        }
        Node<ITEM> *current = this->first;
        Node<ITEM> *otherCurrent = other.first;
        while (current != nullptr) {
            if (current->item != otherCurrent->item) {
                return false;
            }
            current = current->next;
            otherCurrent = otherCurrent->next;
        }
        return true;
    }
    bool operator!=(const Bag<ITEM> &other) const {
        return !(*this == other);
    }
    bool operator<(const Bag<ITEM> &other) const {
        if (this->n < other.n) {
            return true;
        }
        if (this->n > other.n) {
            return false;
        }
        Node<ITEM> *current = this->first;
        Node<ITEM> *otherCurrent = other.first;
        while (current != nullptr) {
            if (current->item < otherCurrent->item) {
                return true;
            }
            if (current->item > otherCurrent->item) {
                return false;
            }
            current = current->next;
            otherCurrent = otherCurrent->next;
        }
        return false;
    }
    bool operator>(const Bag<ITEM> &other) const {
        return other < *this;
    }
    bool operator<=(const Bag<ITEM> &other) const {
        return !(other < *this);
    }
    bool operator>=(const Bag<ITEM> &other) const {
        return !(*this < other);
    }

    ITEM get(int i) {
        Node<ITEM> *current = first;
        for (int j = 0; j < i; j++) {
            current = current->next;
        }
        return current->item;
    }
    ITEM pickRandom() {
        if (n == 0) {
            throw std::out_of_range("Empty bag");
        }
        srand(time(NULL));
        int randomIndex = rand() % n;
        Node<ITEM> *current = first;
        for (int i = 0; i < randomIndex; i++) {
            current = current->next;
        }
        return current->item;
    }
    ITEM removeRandom() {
        if (n == 0) {
            throw std::out_of_range("Empty bag");
        }
        srand(time(NULL));
        int randomIndex = rand() % n;
        Node<ITEM> *current = first;
        for (int i = 0; i < randomIndex; i++) {
            current = current->next;
        }
        ITEM item = current->item;
        if (current == first) {
            first = first->next;
        } else {
            Node<ITEM> *previous = first;
            while (previous->next != current) {
                previous = previous->next;
            }
            previous->next = current->next;
        }
        delete current;
        n--;
        return item;
    }
};
#endif //GRAPH_CPP_BAG_H
