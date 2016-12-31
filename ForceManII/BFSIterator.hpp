#pragma once


/** \brief an iterator that traverses the graph via a breadth first search
 *
 *  The basic idea is on some iteration we are looking at a node, call it i.
 *  We then want to visit all nodes that are connected to i, possibly storing
 *  the paths.
 *
 */
template<typename T>
class BFSIterator{
    const T& graph_;
    size_t max_depth_;
    std::dequeue<std::vector<size_t>> to_go_;
    std::dequeue<std::vector<size_t>> current_;
    void next();
public:
    BFSIterator(const T& graph,size_t start=0):
        graph_(graph),found_({{start}}){
    }

};

void BFSIterator<T>::next(){

}
