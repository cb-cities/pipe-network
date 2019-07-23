//
// Created by Renjie Wu on 2019-07-19.
//

#ifndef PIPE_NETWORK_PUMP_H
#define PIPE_NETWORK_PUMP_H





class Pipe : public Pump {
public:
    //! Constructor with two end nodes, length, diameter roughness and pipe status
    //! \param[in] pipe_prop struct with properties for the pipe
    Pump(const Pipe_prop& pipe_prop)
            : Link(pipe_prop.id, pipe_prop.node1, pipe_prop.node2) {
        pipe_info_["type"] = PIPE;
        pipe_info_["length"] = pipe_prop.length;
        pipe_info_["diameter"] = pipe_prop.diameter;
        pipe_info_["roughness"] = pipe_prop.roughness;
        pipe_info_["status"] = pipe_prop.status;
    };

    //! Virtual destructor
    ~Pipe() override{};

    //! Return link info
    std::map<std::string, double> link_info() const override {
        return pipe_info_;
    }

private:
    // node information, has key : type, elevation, demand, leak_area
    std::map<std::string, double> pipe_info_;
};

#endif //PIPE_NETWORK_PUMP_H
