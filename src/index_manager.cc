#include "index_manager.h"

//! Constructor with an inital value for index
pipenetwork::IndexManager::IndexManager(Index idx) : index_{idx}, next_{idx} {};

//! Current index
pipenetwork::Index pipenetwork::IndexManager::current_index() const {
  return index_;
}

//! Maximum index value
pipenetwork::Index pipenetwork::IndexManager::max() const {
  return std::numeric_limits<Index>::max();
}

//! Generate index
pipenetwork::Index pipenetwork::IndexManager::create_index() {
  std::lock_guard<std::mutex> guard(index_mutex_);
  // Get the next value of index
  // (at init both index and next will be the same)
  index_ = next_;
  // Increment next index by one
  next_ += 1;
  // Return current index
  return index_;
}