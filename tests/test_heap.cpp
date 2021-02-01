#include "id_queue.hpp"
#include "util.hpp"
#include <algorithm>
#include <iostream>

MaxIDQueue<double>
create_small_heap()
{
  MaxIDQueue<double> heap(10);

  heap.push({ 0, 5.0 });
  heap.push({ 1, 6.0 });
  heap.push({ 6, 3.0 });
  heap.push({ 3, 1.0 });
  heap.push({ 4, 9.0 });
  heap.push({ 5, 2.0 });
  heap.push({ 2, 3.0 });
  heap.push({ 7, 8.0 });
  heap.push({ 8, 1.0 });
  heap.push({ 9, 0.0 });

  return heap;
}

void
run_small_heap_test()
{
  auto heap = create_small_heap();
  heap.remove_id(5);

  std::vector<IDKeyPair<double>> removed_points;
  while (!heap.empty()) {
    std::cout << "ID: " << heap.peek().value().id << " Key: " << heap.peek().value().key << std::endl;
    removed_points.push_back({ heap.peek().value().id, heap.peek().value().key });
    heap.pop();
  }

  for (size_t i = 0; i < removed_points.size(); i++)
    heap.push(removed_points[i]);

  std::cout << "Again printing the heap of size " << heap.size() << ":" << std::endl;
  while (!heap.empty()) {
    std::cout << "ID: " << heap.peek().value().id << " Key: " << heap.peek().value().key << std::endl;
    removed_points.push_back({ heap.peek().value().id, heap.peek().value().key });
    heap.pop();
  }
}

void
run_negative_double_values_test(const size_t N, const size_t N_rand_id)
{
  // Generate N negative values and insert them to the heap
  auto heap = MaxIDQueue<double>(N);
  std::vector<double> values_added_to_heap(N);

  for (size_t i = 0; i < N; i++) {
    auto rand_val = Util::rand_double(-100.0, -1.0);
    heap.push({ i, rand_val });
    values_added_to_heap[i] = rand_val;
  }

  // Remove N_rand_id IDs from the heap
  std::vector<size_t> random_IDs_removed_from_heap(N_rand_id);
  for (size_t i = 0; i < N_rand_id; i++) {
    auto rand_id= Util::rand_size_t(0, N-1); // The heap has N-i elements on i^th iteration
    heap.remove_id(rand_id);
    values_added_to_heap[rand_id] = -std::numeric_limits<double>::max();
    random_IDs_removed_from_heap[i] = rand_id;
  }

  std::sort(values_added_to_heap.begin(), values_added_to_heap.end(), std::greater<>());

  std::vector<double> values_retrieved_from_heap;
  while (!heap.empty()) {
    values_retrieved_from_heap.push_back(heap.peek().value().key);
    heap.pop();
  }

  // Compare only the first heap.size() elements of the two arrays
  if (std::equal(values_added_to_heap.begin(), values_added_to_heap.begin() + heap.size(), values_retrieved_from_heap.begin())) {
    std::cout << "OK.";
  } else {
    std::cout << "Not Equal." << std::endl;
    std::cout << "Values from heap: ";
    for (auto i = 0; i < values_retrieved_from_heap.size(); i++) {
      if (values_added_to_heap[i] != values_retrieved_from_heap[i]) {
        std::cout << values_added_to_heap[i] << " " << values_retrieved_from_heap[i] << std::endl;
      }
    }
    std::cout << std::endl;
  }

  // for (auto i: heap_values) {
  //  std::cout << i << " ";
  //}
}

int
main()
{
  for (auto i = 0; i < 100000; i++) {
    run_negative_double_values_test(100, 10);
  }
  return 0;
}