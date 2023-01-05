#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include <stdio.h>
#include <math.h>                   //数学　ヘッダーファイルの読み込み
#include <iostream>                 //出力　ヘッダーファイルの読み込み
#include <iomanip>
#include <vector>
#include <tuple>
#include <map>

#include "vtsimc.hpp"

using namespace::std;
namespace py = pybind11;

PYBIND11_MODULE(vtsimc, m){
    m.attr("SOLVE_LU")   = SOLVE_LU;
    m.attr("SOLVE_SOR")  = SOLVE_SOR;
    
    m.attr("SN_NONE")    = SN_NONE;
    m.attr("SN_CALC")    = SN_CALC;
    m.attr("SN_FIX")     = SN_FIX;
    m.attr("SN_DLY")     = SN_DLY;
    
    m.attr("VN_SIMPLE")  = VN_SIMPLE;
    m.attr("VN_GAP")     = VN_GAP;
    m.attr("VN_FIX")     = VN_FIX;
    m.attr("VN_AIRCON")  = VN_AIRCON;
    m.attr("VN_FAN")     = VN_FAN;
    
    m.attr("TN_SIMPLE")  = TN_SIMPLE;
    m.attr("TN_AIRCON")  = TN_AIRCON;
    m.attr("TN_SOLAR")   = TN_SOLAR;
    m.attr("TN_GROUND")  = TN_GROUND;
    m.attr("TN_HEATER")  = TN_HEATER;
    
    m.attr("AC_NONE")    = AC_NONE;
    m.attr("AC_RAC")     = AC_RAC;
    m.attr("AC_DUCT_C")  = AC_DUCT_C;
    m.attr("AC_RAC2")    = AC_RAC2;
    m.attr("AC_CRIEPI")  = AC_CRIEPI;

    m.attr("AC_AUTO")    = AC_AUTO;
    m.attr("AC_HEATING") = AC_HEATING;
    m.attr("AC_COOLING") = AC_COOLING;
    m.attr("AC_STOP")    = AC_STOP;
    
    py::class_<CalcStatus>(m, "CalcStatus")
        .def(py::init<>())
        .def_readwrite("length",    &CalcStatus::length)
        .def_readwrite("t_step",    &CalcStatus::t_step)
        .def_readwrite("solve",     &CalcStatus::solve)
        .def_readwrite("step_p",    &CalcStatus::step_p) 
        .def_readwrite("vent_err",  &CalcStatus::vent_err) 
        .def_readwrite("step_t",    &CalcStatus::step_t) 
        .def_readwrite("thrm_err",  &CalcStatus::thrm_err)
        .def_readwrite("conv_err",  &CalcStatus::conv_err) 
        .def_readwrite("sor_ratio", &CalcStatus::sor_ratio) 
        .def_readwrite("sor_err",   &CalcStatus::sor_err);
    
    py::class_<VTSim>(m, "VTSim")
        .def(py::init<>())
        .def("setup",              &VTSim::setup,         "")
        .def("sn_add",             &VTSim::sn_add,        "")
        .def("set_node",           &VTSim::set_node,      "")
        .def("vn_add",             &VTSim::vn_add,        "")
        .def("tn_add",             &VTSim::tn_add,        "")
        .def("vn_aircon_add",      &VTSim::vn_aircon_add, "")
        .def("tn_aircon_add",      &VTSim::tn_aircon_add, "")
        .def("calc",               &VTSim::calc,          "")
        .def("result",             &VTSim::result,        "")
        .def_readwrite("sts",      &VTSim::sts,           "")
        .def_readwrite("sn",       &VTSim::sn,            "")
        .def_readwrite("node",     &VTSim::node,          "")
        .def_readwrite("v_net",    &VTSim::v_net,         "")
        .def_readwrite("t_net",    &VTSim::t_net,         "")
        .def_readwrite("vn",       &VTSim::vn,            "")
        .def_readwrite("tn",       &VTSim::tn,            "")
        .def_readwrite("v_idc",    &VTSim::v_idc,         "")
        .def_readwrite("c_idc",    &VTSim::c_idc,         "")
        .def_readwrite("t_idc",    &VTSim::t_idc,         "")
        .def_readwrite("x_idc",    &VTSim::x_idc,         "")
        .def_readwrite("i_vn_ac",  &VTSim::i_vn_ac,       "")
        .def_readwrite("i_tn_ac",  &VTSim::i_tn_ac,       "");

    py::class_<Node>(m, "Node")
        .def(py::init<long, int, tuple<int, int, int, int>>())
        .def_readwrite("i",     &Node::i,     "")
        .def_readwrite("s_i",   &Node::s_i,   "")
        .def_readwrite("flag",  &Node::flag,  "")
        .def_readwrite("p",     &Node::p,     "")
        .def_readwrite("c",     &Node::c,     "")
        .def_readwrite("t",     &Node::t,     "")
        .def_readwrite("x",     &Node::x,     "")
        .def_readwrite("m",     &Node::m,     "")
        .def_readwrite("mx",    &Node::mx,    "")
        .def_readwrite("h_sr",  &Node::h_sr,  "")
        .def_readwrite("h_inp", &Node::h_inp, "")
        .def_readwrite("v",     &Node::v,     "")
        .def_readwrite("beta",  &Node::beta,  "");

    py::class_<Vent_Net>(m, "Vent_Net")
        .def(py::init<long, int, int, int, int, vector<double>, vector<double>>())
        .def_readwrite("i",       &Vent_Net::i,       "")
        .def_readwrite("i1",      &Vent_Net::i1,      "")
        .def_readwrite("i2",      &Vent_Net::i2,      "")
        .def_readwrite("vn_type", &Vent_Net::vn_type, "")
        .def_readwrite("h1",      &Vent_Net::h1,      "")
        .def_readwrite("h2",      &Vent_Net::h2,      "")
        .def_readwrite("alpha",   &Vent_Net::alpha,   "")
        .def_readwrite("area",    &Vent_Net::area,    "")
        .def_readwrite("a",       &Vent_Net::a,       "")
        .def_readwrite("n",       &Vent_Net::n,       "")
        .def_readwrite("qv",      &Vent_Net::qv,      "")
        .def_readwrite("qt",      &Vent_Net::qt,      "")
        .def_readwrite("eta",     &Vent_Net::eta,     "")
        .def_readwrite("q_max",   &Vent_Net::q_max,   "")
        .def_readwrite("p_max",   &Vent_Net::p_max,   "")
        .def_readwrite("q1",      &Vent_Net::q1,      "")
        .def_readwrite("p1",      &Vent_Net::p1,      "");

    py::class_<Thrm_Net>(m, "Thrm_Net")
        .def(py::init<long, int, int, int, int>())
        .def_readwrite("i",              &Thrm_Net::i,              "")
        .def_readwrite("i1",             &Thrm_Net::i1,             "")
        .def_readwrite("i2",             &Thrm_Net::i2,             "")
        .def_readwrite("tn_type",        &Thrm_Net::tn_type,        "")
        .def_readwrite("cdtc",           &Thrm_Net::cdtc,           "")
        .def_readwrite("ms",             &Thrm_Net::ms,             "")
        .def_readwrite("area",           &Thrm_Net::area,           "")
        .def_readwrite("phi_0",          &Thrm_Net::phi_0,          "")
        .def_readwrite("cof_r",          &Thrm_Net::cof_r,          "")
        .def_readwrite("cof_phi",        &Thrm_Net::cof_phi,        "")
        .def_readwrite("t_dash_gs",      &Thrm_Net::t_dash_gs,      "")
        .def_readwrite("qt",             &Thrm_Net::qt,             "")
        .def_readwrite("aircon_on",      &Thrm_Net::aircon_on,      "")
        .def_readwrite("ac_model",       &Thrm_Net::ac_model,       "")
        .def_readwrite("ac_mode",        &Thrm_Net::ac_mode,        "")
        .def_readwrite("pre_tmp",        &Thrm_Net::pre_tmp,        "")
        .def_readwrite("pre_x",          &Thrm_Net::pre_x,          "")
        .def_readwrite("q_rtd_C",        &Thrm_Net::q_rtd_C,        "")
        .def_readwrite("q_max_C",        &Thrm_Net::q_max_C,        "")
        .def_readwrite("e_rtd_C",        &Thrm_Net::e_rtd_C,        "")
        .def_readwrite("q_rtd_H",        &Thrm_Net::q_rtd_H,        "")
        .def_readwrite("q_max_H",        &Thrm_Net::q_max_H,        "")
        .def_readwrite("e_rtd_H",        &Thrm_Net::e_rtd_H,        "")
        .def_readwrite("dualcompressor", &Thrm_Net::dualcompressor, "")
        .def_readwrite("q_hs_rtd_H",     &Thrm_Net::q_hs_rtd_H,     "")
        .def_readwrite("q_hs_mid_H",     &Thrm_Net::q_hs_mid_H,     "")
        .def_readwrite("q_hs_min_H",     &Thrm_Net::q_hs_min_H,     "") 
        .def_readwrite("P_hs_rtd_H",     &Thrm_Net::P_hs_rtd_H,     "") 
        .def_readwrite("V_fan_rtd_H",    &Thrm_Net::V_fan_rtd_H,    "") 
        .def_readwrite("V_fan_mid_H",    &Thrm_Net::V_fan_mid_H,    "") 
        .def_readwrite("P_fan_rtd_H",    &Thrm_Net::P_fan_rtd_H,    "") 
        .def_readwrite("P_fan_mid_H",    &Thrm_Net::P_fan_mid_H,    "") 
        .def_readwrite("V_hs_dsgn_H",    &Thrm_Net::V_hs_dsgn_H,    "")
        .def_readwrite("q_hs_rtd_C",     &Thrm_Net::q_hs_rtd_C,     "")  
        .def_readwrite("q_hs_mid_C",     &Thrm_Net::q_hs_mid_C,     "") 
        .def_readwrite("q_hs_min_C",     &Thrm_Net::q_hs_min_C,     "")  
        .def_readwrite("P_hs_rtd_C",     &Thrm_Net::P_hs_rtd_C,     "")
        .def_readwrite("V_fan_rtd_C",    &Thrm_Net::V_fan_rtd_C,    "")
        .def_readwrite("V_fan_mid_C",    &Thrm_Net::V_fan_mid_C,    "")
        .def_readwrite("P_fan_rtd_C",    &Thrm_Net::P_fan_rtd_C,    "")
        .def_readwrite("P_fan_mid_C",    &Thrm_Net::P_fan_mid_C,    "")
        .def_readwrite("V_hs_dsgn_C",    &Thrm_Net::V_hs_dsgn_C,    "")
        .def_readwrite("To",             &Thrm_Net::To,             "")
        .def_readwrite("ho",             &Thrm_Net::ho,             "");
}          