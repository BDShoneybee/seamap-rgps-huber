Dear Editors of Satellite Navigation,

Please consider our original research manuscript, “Causal robust trend preprocessing of marine relative-GPS range and bearing: reproducible evidence and operational limits,” by Menghui Fan, Shouchuan Fang and Fujiang Liu.

The study addresses a scalar observation interface in marine relative-GPS navigation. It reconstructs an engineering Huber trend preprocessor from source code, states its causal initialization and missing-data rules, and evaluates its gains and limitations through reproducible Monte Carlo experiments, polar geometry, component ablations and chronological replay of ship-collected Seamap observations. A portable C++ implementation and Python harness are numerically cross-checked.

The manuscript explicitly separates simulated truth-based errors from real-data continuity and prediction diagnostics. We do not claim measured absolute accuracy in the absence of an independent trajectory. Negative bearing-tracking results and abrupt-change lag are reported alongside synthetic outlier improvements. We believe this transparent treatment of robust processing at a marine navigation interface is relevant to the journal's navigation signal processing and industrial applications readership.

Data, code, results and LaTeX sources: [seamap-rgps-huber](https://github.com/BDShoneybee/seamap-rgps-huber).

The authors declare no competing interests. Funding is currently left blank at the authors' request.

**Author completion before submission:** insert confirmation that all authors approve the final manuscript and that its content has not been published or submitted elsewhere. Confirm funding or explicitly declare no funding if appropriate. These statements have not been assumed by the drafting assistant. The review copy removes authors and the identifying repository URL; supporting material should be supplied in accordance with the journal's double-anonymous review requirements.

Sincerely,

Shouchuan Fang, corresponding author  
BGP Offshore, BGP Inc., CNPC, China  
fangshouchuan@cnpc.com.cn
