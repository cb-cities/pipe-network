#ifndef PIPE_NETWORK_MESH_H_
#define PIPE_NETWORK_MESH_H_

#include <cmath>

#include <array>
#include <exception>
#include <iostream>
#include <map>
#include <memory>
#include <vector>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/betweenness_centrality.hpp>
#include <boost/graph/graph_traits.hpp>

#include "node.h"
#include "pipe.h"

//! Mesh class
//! \brief Class for mesh that contains node and pipe pointers
class Mesh {

  typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS>
      graph;
  typedef boost::graph_traits<graph>::vertex_descriptor vertex;
  typedef boost::graph_traits<graph>::edge_descriptor edge;
  typedef boost::property_map<graph, boost::vertex_index_t>::type vindexmap;
  typedef std::map<edge, int> stdeindex;
  typedef boost::associative_property_map<stdeindex> eindexmap;

 public:
  // Constructor with id
  //! \param[in] id mesh id
  Mesh(unsigned id) : id_{id} {};

  //! Destructor
  ~Mesh() {
    nodes_.clear();
    pipes_.clear();
  };

  //! Return id
  //! \retval id_ id of the mesh
  unsigned id() const { return id_; }

  //! Create a nodal pointer based on input id and coordinates
  //! \param[in] id the node id
  //! \param[in] coords the coordinates of the node
  void create_node(unsigned id, Eigen::Vector3d coords) {
    nodes_.emplace_back(std::make_shared<pipenetwork::Node>(id, coords));
  }

  //! Check whether there are nodes with duplicate indices or coordinates
  void duplicate_node() {
    for (unsigned i = 0; i < nodes_.size() - 1; ++i) {
      for (unsigned j = i + 1; j < nodes_.size(); ++j) {
        if (nodes_.at(i)->id() == nodes_.at(j)->id()) {
          std::cout << "node id: " << nodes_.at(i)->id() << '\n';
          throw std::runtime_error(
              "duplicate nodal indices found, check input");
        }
        if (nodes_.at(i)->coordinates() == nodes_.at(j)->coordinates()) {
          std::cout << "node id: " << nodes_.at(i)->id() << " and "
                    << nodes_.at(j)->id() << '\n';
          throw std::runtime_error(
              "duplicate nodal coordinates found, check input");
        }
      }
    }
  }

  //! Return a nodal pointer for a given nodal index
  //! \param[in] id the nodal index
  //! \retval nodeptr the nodal pointer with the given nodal index
  std::shared_ptr<pipenetwork::Node> node(unsigned id) {
    std::shared_ptr<pipenetwork::Node> nodeptr = nullptr;
    for (const auto& node : nodes_) {
      if (node->id() == id) {
        nodeptr = node;
        break;
      }
    }
    if (nodeptr == nullptr) {
      std::cout << "node not found, id: " << id << '\n';
      throw std::runtime_error(
          "node with indicates id is not found, check input");
    } else
      return nodeptr;
  }

  //! Create a pipe pointer based on input idices of the pipe and the nodes at
  //! its ends 
  //! \param[in] pipeid index of the pipe 
  //! \param[in] nodeid1 and nodeid2 indices of the nodes at pipe ends
  void create_pipe(unsigned pipeid, unsigned nodeid1, unsigned nodeid2) {
    std::array<std::shared_ptr<pipenetwork::Node>, 2> nodes;
    nodes.at(0) = node(nodeid1);
    nodes.at(1) = node(nodeid2);
    nodes_of_pipe_.insert(
        std::pair<unsigned, std::array<std::shared_ptr<pipenetwork::Node>, 2>>(
            pipeid, nodes));
    pipes_.emplace_back(std::make_unique<pipenetwork::Pipe>(pipeid, nodes));
  }

  //! Check whether there are pipes with duplicate indices
  void duplicate_pipe() {
    for (unsigned i = 0; i < pipes_.size() - 1; ++i) {
      for (unsigned j = i + 1; j < pipes_.size(); ++j) {
        if (pipes_.at(i)->id() == pipes_.at(j)->id()) {
          std::cout << "pipe id: " << pipes_.at(i)->id() << '\n';
          throw std::runtime_error("duplicate pipe indices found, check input");
        }
      }
    }
  }

  //! Return a pipe pointer for a given pipe index
  //! \param[in] id the pipe index
  //! \retval pipeptr the pipe pointer with the given pipe index
  //  std::unique_ptr<pipenetwork::Pipe> pipe(unsigned id) {
  //	  std::unique_ptr<pipenetwork::Pipe> pipeptr = nullptr;
  //	  for(const auto& pipe : pipes_) {
  //		  if(pipe->id() == id) {
  //			  pipeptr = pipe;
  //			  break;
  //		  }
  //	  }
  //	  if(pipeptr == nullptr)
  //		  throw std::runtime_error("pipe with indicates id is not found,
  //check input"); 	  else 		  return pipeptr;
  //}

  //! Return the number of nodes in the mesh
  //! \retval nodes_.size() number of nodes
  unsigned nnodes() const { return nodes_.size(); }

  //! Return the number of pipes in the mesh
  //! \retval pipes_.size() number of pipes
  unsigned npipes() const { return pipes_.size(); }

  //! Check whether redundant node exists
  void redundant_node() {
    for (const auto& node : nodes_) {
      bool connect = false;
      for (const auto& pipe : pipes_) {
        if (nodes_of_pipe_.at(pipe->id()).at(0)->id() == node->id() ||
            nodes_of_pipe_.at(pipe->id()).at(1)->id() == node->id()) {
          connect = true;
          break;
        }
      }
      if (connect == false)
        throw std::runtime_error("redundant node exist, check input");
    }
  }

  //! Return coordinates of all the nodes in the mesh
  //! \retval nodal_coordinates coordinates of all the nodes
  std::vector<Eigen::Vector3d> nodal_coordinates() {
    std::vector<Eigen::Vector3d> nodal_coordinates;
    for (const auto& node : nodes_)
      nodal_coordinates.emplace_back(node->coordinates());
    return nodal_coordinates;
  }

  //! Return degree centrality (number of pipe connected to the node) of a given
  //! node param[in] id index of the interested node 
  //! \retval degree_centrality degree centrality of the node
  unsigned degree_centrality(unsigned id) {
    unsigned degree_centrality = 0;
    for (const auto& pipe : pipes_)
      if (id == nodes_of_pipe_.at(pipe->id()).at(0)->id() ||
          id == nodes_of_pipe_.at(pipe->id()).at(1)->id())
        degree_centrality++;
    return degree_centrality;
  }

  //! Calculate the betweenness centrality of each node and pipe
  void calc_betweenness_centrality() {
    //	  typedef boost::adjacency_list<boost::vecS, boost::vecS,
    //boost::undirectedS> graph;
    graph g;

    //	  typedef boost::graph_traits<graph>::vertex_descriptor vertex;
    //	  std::map<int, vertex> vertexid;
    std::vector<vertex> vertex_vec;
    for (unsigned i = 0; i < nodes_.size(); ++i) {
      vertex_vec.emplace_back(boost::add_vertex(g));
      vertexid_.insert(
          std::pair<int, vertex>(nodes_.at(i)->id(), vertex_vec.at(i)));
    }

    //	  typedef boost::graph_traits<graph>::edge_descriptor edge;
    //	  std::map<int, edge> edgeid;
    std::vector<edge> edge_vec;
    for (unsigned i = 0; i < pipes_.size(); ++i) {
      edge_vec.emplace_back(
          (boost::add_edge(
               vertexid_.at(nodes_of_pipe_.at(pipes_.at(i)->id()).at(0)->id()),
               vertexid_.at(nodes_of_pipe_.at(pipes_.at(i)->id()).at(1)->id()),
               g))
              .first);
      edgeid_.insert(std::pair<int, edge>(pipes_.at(i)->id(), edge_vec.at(i)));
    }

    //	  typedef boost::property_map<graph, boost::vertex_index_t>::type
    //vindexmap;
    vindexmap v_index = get(boost::vertex_index, g);
    std::vector<double> v_centrality_vector(boost::num_vertices(g), 0.0);
    boost::iterator_property_map<std::vector<double>::iterator, vindexmap>
        v_centrality_map(v_centrality_vector.begin(), v_index);

    //	  typedef std::map<edge, int> stdeindex;
    stdeindex std_e_index;
    //	  typedef boost::associative_property_map<stdeindex> eindexmap;
    eindexmap e_index(std_e_index);
    for (int i = 0; i < edge_vec.size(); ++i)
      std_e_index.insert(std::pair<edge, int>(edge_vec.at(i), i));
    std::vector<double> e_centrality_vector(boost::num_edges(g), 0.0);
    boost::iterator_property_map<std::vector<double>::iterator, eindexmap>
        e_centrality_map(e_centrality_vector.begin(), e_index);

    boost::brandes_betweenness_centrality(g, v_centrality_map,
                                          e_centrality_map);
    node_centrality_map_ = v_centrality_map;
    pipe_centrality_map_ = e_centrality_map;
    is_betweenness_centrality_ = true;
  }

  //! Retrun betweenness centrality of indicated node
  //! \param[in] id index of the interested node
  //! \retval betweenness centrality of the node
  double betweenness_centrality_node(unsigned id) {
    if (is_betweenness_centrality_ = false)
      throw std::runtime_error("the value hasn't been calculated");
    return node_centrality_map_[vertexid_.at(id)];
  }

  //! Retrun betweenness centrality of indicated pipe
  //! \param[in] id index of the interested pipe
  //! \retval betweenness centrality of the pipe
  double betweenness_centrality_pipe(unsigned id) {
    if (is_betweenness_centrality_ = false)
      throw std::runtime_error("the value hasn't been calculated");
    return pipe_centrality_map_[edgeid_.at(id)];
  }

 private:
  //! mesh id
  unsigned id_{std::numeric_limits<unsigned>::max()};
  //! nodal pointers
  std::vector<std::shared_ptr<pipenetwork::Node>> nodes_;
  //! pipe pointers
  std::vector<std::unique_ptr<pipenetwork::Pipe>> pipes_;
  //! corelate pipe id with nodes at its ends
  std::map<unsigned, std::array<std::shared_ptr<pipenetwork::Node>, 2>>
      nodes_of_pipe_;
  //! calculation status of betweenness centrality
  bool is_betweenness_centrality_{false};
  //! map between node id and vertex iterator
  std::map<int, vertex> vertexid_;
  //! map between pipe id and edge iterator
  std::map<int, edge> edgeid_;
  //! node centrality map
  boost::iterator_property_map<std::vector<double>::iterator, vindexmap>
      node_centrality_map_;
  //! pipe centrality map
  boost::iterator_property_map<std::vector<double>::iterator, eindexmap>
      pipe_centrality_map_;
};

#endif  // PIPE_NETWORK_MESH_H_
