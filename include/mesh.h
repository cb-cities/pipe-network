#ifndef PIPE_NETWORK_MESH_H
#define PIPE_NETWORK_MESH_H
#include <fstream>
#include <iomanip>

#include "mesh_components.h"
namespace pipenetwork {

////! Mesh class
////! \brief Class for mesh that contains node and pipe pointers
class Mesh {
 public:
  //! Constructor with id
  //! \param[in] id mesh id
  explicit Mesh(std::string& name) : name_{name} {};

  //! Return name
  //! \retval name_ name of the mesh
  std::string name() const { return name_; }

  //  //! Create mesh from input object
  //  //! \param[in] IO pointer to the input object
  //  void create_mesh_from_io(const std::shared_ptr<IO>& IO);

  //! Create junction pointers
  //! \param[in] junc_props vector of junction properties
  void create_nodes(const std::vector<JunctionProp>& junc_props,
                    const std::vector<ReservoirProp>& res_props) {
    mesh_nodes_ = std::make_shared<MeshNodes>(junc_props, res_props);
    find_leak_nids();
  };

  //! Create Reservoir pointers
  //! \param[in] res_props vector of reservoir properties
  void create_links(const std::vector<PipeProp>& pipe_props,
                    const std::vector<PumpProp>& pump_props,
                    const std::vector<ValveProp>& valve_props) {
    mesh_links_ = std::make_shared<MeshLinks>(pipe_props, pump_props,
                                              valve_props, *mesh_nodes_);
  };

  //! Create graph for mesh
  void create_mesh_graph() {
    mesh_graph_ = std::make_shared<MeshGraph>(mesh_nodes_, mesh_links_);
    find_iso_components_();
  }

  //! get nodes
  std::shared_ptr<MeshNodes> nodes() const { return mesh_nodes_; }
  //! get links
  std::shared_ptr<MeshLinks> links() const { return mesh_links_; }
  //! get links
  std::shared_ptr<MeshGraph> mesh_graph() const { return mesh_graph_; }
  //! get isolated nodes
  const std::vector<Index>& iso_nodes() const { return iso_nodes_; }
  //! get isolated links
  const std::vector<Index>& iso_links() const { return iso_links_; }
  //! get leak nids
  const std::vector<Index>& leak_nids() const { return leak_nids_; }

  //! Print summary for the mesh
  void print_summary();

 private:
  //! mesh name
  std::string name_;
  //! mesh nodes
  std::shared_ptr<MeshNodes> mesh_nodes_;
  //! mesh links
  std::shared_ptr<MeshLinks> mesh_links_;
  //! graph
  std::shared_ptr<MeshGraph> mesh_graph_;
  //! isolated nodes ids
  std::vector<Index> iso_nodes_;
  //! isolated link ids
  std::vector<Index> iso_links_;
  //! leak nids
  std::vector<Index> leak_nids_;
  //! Find isolated junctions and links for the mesh
  void find_iso_components_();
  //! find isolated nodes
  void find_iso_nodes_();
  //! find isolated links
  void find_iso_links_();
  //! find leak node lids
  void find_leak_nids();
};

}  // namespace pipenetwork

#endif  // PIPE_NETWORK_MESH_H
