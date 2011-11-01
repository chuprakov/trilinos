#ifndef MUELU_LINKEDLIST_HPP
#define MUELU_LINKEDLIST_HPP

#include "MueLu_ConfigDefs.hpp"

/* ------------------------------------------------------------------------- */
/* linked list structures from ML for holding free node information          */
/* ------------------------------------------------------------------------- */

namespace MueLu {

  typedef struct MueLu_Node_Struct
  {
    int nodeId;
    struct MueLu_Node_Struct *next;
  } MueLu_Node;
  
  class LinkedList {

  public:    
    LinkedList() : nodeHead(NULL), nodeTail(NULL) { }

    ~LinkedList() {
      while (nodeHead != NULL)
        DeleteHead();
    }

    bool IsEmpty() {
      return nodeHead == NULL;
    }

    void Add(int iNode) {
      MueLu_Node *newNode = new MueLu_Node;
      newNode->nodeId = iNode;
      newNode->next = NULL;
      if (nodeHead == NULL) {
          nodeHead = newNode;
          nodeTail = newNode;
        } else {
          nodeTail->next = newNode;
          nodeTail = newNode;
        }
    }

    int Pop() { // get head and remove first node
      if (IsEmpty()) throw(1);

      int iNode = nodeHead->nodeId;
      DeleteHead();
      return iNode;
    }

  private:
    MueLu_Node *nodeHead;
    MueLu_Node *nodeTail;

    void DeleteHead() {
      if (IsEmpty()) throw(1);
      
      MueLu_Node *newNode = nodeHead;
      nodeHead = newNode->next;
      delete newNode;
    }

  };

}
#endif
//TODO: nodeTail unused -> remove?
