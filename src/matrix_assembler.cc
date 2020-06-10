#include "matrix_assembler.h"

// pipenetwork::MatrixAssembler::MatrixAssembler(
//        const std::shared_ptr<pipenetwork::Mesh>& mesh,
//        std::shared_ptr<pipenetwork::Curves>& curves_info, bool pdd_mode)
//        : mesh_{mesh}, curves_info_{curves_info}, pdd_{pdd_mode} {
//
//    init_variable_vectors();
////    assemble_balance_headloss_matrix();
////    initialize_jacobian();
//
//}

// Initialize variable vector
// 0 to nnode-1 element: nodal head
// nnode to 2*nnode-1 element: nodal demand
// 2*nnode to 2*nnode+nlinks-1 element: pipe discharge
// 2*nnode+nlinks to 2*nnode+nlinks+nleaks element: leak discharge
void pipenetwork::linear_system::Variables::init_variable_vectors() {
  auto nodes = mesh_->nodes();
  auto links = mesh_->links();
  auto leak_ids = mesh_->leak_nids();
  Index nnodes = nodes->nnodes();
  Index nlinks = links->nlinks();
  Index nleaks = leak_ids.size();

  resize_init_variables(nnodes, nlinks, nleaks);

  init_nodes_vecs(nodes);
  init_link_vecs(links, 2 * nnodes);
  init_leak_nodes(nodes, leak_ids, 2 * nnodes + nlinks);
}

void pipenetwork::linear_system::Variables::resize_init_variables(
    Index nnodes, Index nlinks, Index nleaks) {

  variable_vec_.resize(2 * nnodes + nlinks + nleaks);
  demands_heads_vec_.resize(nnodes);
  elevations_.resize(nnodes);
  link_resistance_coeff_vec_.resize(nlinks);
  link_minor_loss_coeff_vec_.resize(nlinks);
  leak_areas_.resize(nleaks);
}

void pipenetwork::linear_system::Variables::init_nodes_vecs(
    const std::shared_ptr<pipenetwork::MeshNodes>& nodes) {
  auto junctions_map = nodes->junctions();
  auto reservors_map = nodes->reservoirs();
  auto nnodes = nodes->nnodes();

  for (const auto& index_junction : junctions_map) {
    auto nid = index_junction.first;
    auto junction = index_junction.second;
    auto junc_prop = junction->property();
    demands_heads_vec_[nid] = junc_prop.demand;
    variable_vec_[nnodes + nid] = junc_prop.demand;

    elevations_[nid] = junc_prop.elevation;
    variable_vec_[nid] = junc_prop.elevation;
  }

  for (const auto& index_res : reservors_map) {
    auto nid = index_res.first;
    auto res = index_res.second;
    demands_heads_vec_[nid] = res->head();
    variable_vec_[nid] = res->head();
    elevations_[nid] = res->head();
    variable_vec_[nnodes + nid] = 0;
  }
}

void pipenetwork::linear_system::Variables::init_link_vecs(
    const std::shared_ptr<pipenetwork::MeshLinks>& links,
    Index link_start_idx) {
  auto pipes_map = links->pipes();
  auto pumps_map = links->pumps();
  auto valves_map = links->valves();
  auto links_map = links->links();

  update_minor_resis_coeff(pipes_map);
  update_minor_resis_coeff(pumps_map);
  update_minor_resis_coeff(valves_map);

  for (const auto& index_link : links_map) {
    auto lid = index_link.first;
    variable_vec_[link_start_idx + lid] = INIT_FLOWRATE;
  }
}

template <typename LinkMap>
void pipenetwork::linear_system::Variables::update_minor_resis_coeff(
    const LinkMap& linkmap) {

  for (const auto& index_link : linkmap) {
    auto lid = index_link.first;
    auto link = index_link.second;
    auto link_res_coeff = get_link_res_coeff(link);
    auto minor_loss_coeff = get_link_minor_coeff(link);

    link_resistance_coeff_vec_[lid] = link_res_coeff;
    link_minor_loss_coeff_vec_[lid] = minor_loss_coeff;
  }
}

double pipenetwork::linear_system::Variables::get_link_res_coeff(
    const std::shared_ptr<pipenetwork::Pipe>& pipe) {
  auto pipe_property = pipe->property();
  auto link_res_coeff = HW_COEFF * std::pow(pipe_property.roughness, -1.852) *
                        std::pow(pipe_property.diameter, -4.871) *
                        pipe_property.length;
  return link_res_coeff;
}

double pipenetwork::linear_system::Variables::get_link_res_coeff(
    const std::shared_ptr<pipenetwork::Valve>& valve) {

  auto valve_property = valve->property();
  double link_res_coeff = 0;
  if (valve_property.type == ValveType::TCVALVE) {
    link_res_coeff =
        8 * valve_property.setting /
        (G * std::pow(PI, 2) * std::pow(valve_property.diameter, 4));
  }
  return link_res_coeff;
}

void pipenetwork::linear_system::Variables::init_leak_nodes(
    const std::shared_ptr<pipenetwork::MeshNodes>& nodes,
    const std::vector<Index>& leak_ids, Index leak_start_idx) {
  auto junction_map = nodes->junctions();
  Index i = 0;
  for (const auto& leak_id : leak_ids) {
    auto leak_junction = junction_map.at(leak_id);
    leak_areas_[i] = leak_junction->leak_area();
    variable_vec_[leak_start_idx + i] = INIT_FLOWRATE;
    i++;
  }
}
