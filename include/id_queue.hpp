#pragma once

#include <algorithm>
#include <cassert>
#include <vector>

#include <iostream>
#include <optional>

#include "constants.hpp"

template<typename T>
struct IDKeyPair
{
  size_t id;
  T key;
};

// A priority queue where the elements are IDs from 0 to id_count-1 where id_count is a number that
// is set in the constructor.
template<typename T>
class MaxIDQueue
{
private:
  static constexpr size_t tree_arity = 4;
  static constexpr size_t invalid_id = std::numeric_limits<size_t>::max();

public:
  // MaxIDQueue():heap_size(0){}

  explicit MaxIDQueue(size_t id_count)
    : id_pos(id_count, invalid_id)
    , heap(id_count)
    , heap_size(0)
  {}

  // Returns whether the queue is empty. Equivalent to checking whether size() returns 0.
  bool empty() const { return heap_size == 0; }

  size_t size() const { return heap_size; }

  size_t id_count() const { return id_pos.size(); }

  // Checks whether an element is in the queue.
  bool contains_id(size_t id)
  {
    assert(id < id_count());
    return id_pos[id] != invalid_id;
  }

  // Removes all elements from the queue.
  void clear()
  {
    for (size_t i = 0; i < heap_size; ++i)
      id_pos[heap[i].id] = invalid_id;
    heap_size = 0;
  }

  friend void swap(MaxIDQueue& l, MaxIDQueue& r)
  {
    using std::swap;
    swap(l.id_pos, r.id_pos);
    swap(l.heap, r.heap);
    swap(l.heap_size, r.heap_size);
  }

  // Returns the current key of an element.
  // Undefined if the element is not part of the queue.
  unsigned get_key(size_t id) const
  {
    assert(id < id_count());
    assert(id_pos[id] != invalid_id);
    return heap[id_pos[id]].key;
  }

  // Returns the largest element key pair without removing it from the queue.
  std::optional<IDKeyPair<T>> peek() const
  {
    if (!empty()) {
      return heap.front();
    }
    return std::nullopt;
  }

  // Returns the largest element key pair and removes it form the queue.
  std::optional<IDKeyPair<T>> pop()
  {
    if (!empty()) {
      --heap_size;
      std::swap(heap[0].key, heap[heap_size].key);
      std::swap(heap[0].id, heap[heap_size].id);
      id_pos[heap[0].id] = 0;
      id_pos[heap[heap_size].id] = invalid_id;

      move_down_in_tree(0);
      return heap[heap_size];
    }
    return std::nullopt;
  }

  std::optional<IDKeyPair<T>> remove_id(const size_t id)
  {
    assert(!empty());
    assert(id < id_count());

    if (contains_id(id)) {
      --heap_size;

      size_t pos = id_pos[id];
      std::swap(heap[pos].key, heap[heap_size].key);
      std::swap(heap[pos].id, heap[heap_size].id);
      id_pos[heap[pos].id] = pos;
      id_pos[heap[heap_size].id] = invalid_id;

      size_t parent = (pos - 1) / tree_arity;

      if (pos != 0 && std::isgreater((heap[pos].key - heap[parent].key), DwellRegions::Constants::tolerance)) {
        move_up_in_tree(pos);
      } else {
        move_down_in_tree(pos);
      }

      return heap[heap_size];
    }

    return std::nullopt;
  }

  // Inserts a element key pair.
  // Undefined if the element is part of the queue.
  void push(IDKeyPair<T> p)
  {
    assert(p.id < id_count());
    assert(!contains_id(p.id));

    size_t pos = heap_size;
    ++heap_size;
    heap[pos] = p;
    id_pos[p.id] = pos;
    move_up_in_tree(pos);
  }

  // Updates the key of an element if the new key is smaller than the old key.
  // Does nothing if the new key is larger.
  // Undefined if the element is not part of the queue.
  bool decrease_key(IDKeyPair<T> p)
  {
    assert(p.id < id_count());
    assert(contains_id(p.id));

    size_t pos = id_pos[p.id];

    if (std::isgreater((heap[pos].key - p.key), DwellRegions::Constants::tolerance)) {
      heap[pos].key = p.key;
      move_up_in_tree(pos);
      return true;
    } else {
      return false;
    }
  }

  // Updates the key of an element if the new key is larger than the old key.
  // Does nothing if the new key is smaller.
  // Undefined if the element is not part of the queue.
  bool increase_key(IDKeyPair<T> p)
  {
    assert(p.id < id_count());
    assert(contains_id(p.id));

    size_t pos = id_pos[p.id];

    if (std::isless((heap[pos].key - p.key), -DwellRegions::Constants::tolerance)) {
      heap[pos].key = p.key;
      move_down_in_tree(pos);
      return true;
    } else {
      return false;
    }
  }

private:
  void move_up_in_tree(size_t pos)
  {
    while (pos != 0) {
      size_t parent = (pos - 1) / tree_arity;
      // If the parent is smaller than the key at pos
      if (std::isless((heap[parent].key - heap[pos].key), -DwellRegions::Constants::tolerance)) {
        std::swap(heap[pos], heap[parent]);
        std::swap(id_pos[heap[pos].id], id_pos[heap[parent].id]);
      }
      pos = parent;
    }
  }

  void move_down_in_tree(size_t pos)
  {
    for (;;) {
      size_t first_child = tree_arity * pos + 1;
      if (first_child >= heap_size)
        return; // no children
      size_t largest_child = first_child;
      for (size_t c = first_child + 1; c < std::min(tree_arity * pos + tree_arity + 1, heap_size); ++c) {
        if (std::isless((heap[largest_child].key - heap[c].key), -DwellRegions::Constants::tolerance)) {
          largest_child = c;
        }
      }

      if (std::islessequal((heap[largest_child].key - heap[pos].key), -DwellRegions::Constants::tolerance))
        return; // no child is larger

      std::swap(heap[pos], heap[largest_child]);
      std::swap(id_pos[heap[pos].id], id_pos[heap[largest_child].id]);
      pos = largest_child;
    }
  }

  std::vector<size_t> id_pos;
  std::vector<IDKeyPair<T>> heap;

  size_t heap_size;
};
