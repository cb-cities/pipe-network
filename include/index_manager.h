#ifndef PIPE_NETWORK_INDEX_MANAGER_H
#define PIPE_NETWORK_INDEX_MANAGER_H
#include <limits>
#include <mutex>

#include "settings.h"

namespace pipenetwork {
//! Network index class
//! \brief Base class of a network index
class IndexManager {
 public:
  //! Constructor with an inital value for index
  explicit IndexManager(Index idx = 0);

  //! Current index
  Index current_index() const;

  //! Max value
  Index max() const;

  //! Create index
  Index create_index();

 private:
  //! index
  Index index_{0};
  //! Next value
  Index next_{0};
  //! Mutex
  std::mutex index_mutex_;
};
}  // namespace pipenetwork

#endif  // PIPE_NETWORK_INDEX_MANAGER_H
