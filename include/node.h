#ifndef PIPE_NETWORK_NODE_H_
#define PIPE_NETWORK_NODE_H_

#include <Eigen/Sparse>

//! Node class
//! \brief Class that stores the information about nodes
class Node {

	public:

		// Constructor with id and coordinates
		//! \param[in] id assign as the id_ of the node
		//! \param[in] coords coordinates of the node
		Node(unsigned id, const Eigen::Vector3d& coords) : id_{id} {
			for (unsigned i = 0; i < coords.size(); i++) coordinates_[i] = coords(i);
		}

		//! Destructor
		~Node() = default;

		//! Copy constructor
		Node(const Node&) = delete;

		//! Assignment operator
		Node& operator=(const Node&) = delete;

		//! Move constructor
		Node(Node&&) = delete;

		//! Return id
		//! \param[out] id_ return id of the node
		unsigned id() const { return id_; }

		//! Return coordinates
		//! \param[out] coordinates_ return coordinates of the node
		Eigen::Vector3d coordinates() const { return coordinates_; }

		//! Return coordinates in a particular direction
		//! \param[in] dir direction
		//! \param[out] coordinates_[dir] return coordinates of the node in a particular direction
		double coord_at_dir(unsigned dir) const { return coordinates_[dir]; }

		//! Return number of connection
		//! \param[out] num_of_connection_ return number of connection to the node
		unsigned num_of_connection() const { return num_of_connection_; }

		//! Assign hydraulic head at the node
		//! \param[in] head hydraulic head at the node
		void set_head(double head) {
			head_ = head;
			ishead_ = true;
		}

		//! Return hydraulic head
		//! \param[out] head_ return hydraulic head at the node
		double head() const { return head_; }

		//! Return head assignment status
		//! \param[out] ishead_ return head assignment status at the node
		bool ishead() const {return ishead_; }

		//! Assign discharge at the node
		//! \param[in] discharge discharge at the node
		void set_discharge(double discharge) {
			discharge_ = discharge;
			isdischarge_ = true;
		}

		//! Return discharge
		//! \param[out] discharge_ return discharge at the node
		double discharge() const { return discharge_;}

		//! Return discharge assignment status
		//! \param[out] isdischarge return discharge assignment status at the node
		bool isdischarge() const { return isdischarge_; }

	private:
		//! node id
		unsigned id_{std::numeric_limits<unsigned>::max()};
		//! nodal coordinates
		Eigen::Vector3d coordinates_;
		//! number of connection to the node
		unsigned num_of_connection_{std::numeric_limits<unsigned>::max()};
		//! hydraulic head
		double head_{std::numeric_limits<double>::max()};
		//! discharge
		double discharge_{std::numeric_limits<double>::max()};
		//! whether head is assigned
		bool ishead_{false};
		//! whether discharge is assigned
		bool isdischarge_{false};
};

#endif  // PIPE_NETWORK_NODE_H_
