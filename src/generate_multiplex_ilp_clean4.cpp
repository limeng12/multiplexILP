#include <Rcpp.h>
#include <vector>
#include <unordered_map>
#include <string>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <set>
#include <iomanip>
#include <queue>  // 添加这个头文件
#include <functional>  // 用于greater
#include <string> 
using namespace Rcpp;
using namespace std;

//[[Rcpp::export(rng = false)]]
List generate_correct_multiplex_ilp_locus_C_fixed_binary(
    StringVector primers_5,
    StringVector primers_3,
    StringVector ids_5,
    StringVector ids_3,
    StringVector str_loci_all,
    NumericVector dis_5,
    NumericVector dis_3,
    NumericVector temp_5,
    NumericVector temp_3,
    NumericVector weights,
    NumericVector amplicon_size_low,
    NumericVector amplicon_size_high,
    StringVector compatible_pairs,
    NumericVector num_channels_v,
    NumericVector size_gap_v = 10,
    NumericVector lower_bound = NumericVector::create(0)
  
) {
  
  
  
  
  int num_channels=num_channels_v.at(0);
  double distance_weight = 0.1;
  int size_gap =size_gap_v.at(0);
  double M = 900;
  
  
  
  int n = primers_5.size();
  
  
  if(n>30000){
    Rcpp::Rcout << " too many primers  " << "\n"; 
    
  }
  
  // 1. 建立引物对信息
  struct PrimerPair {
    string pair_id;
    string locus;
    string id5, id3;
    string primers5,primers3;
    string temp5,temp3;
    double weight;
    double distance;
    double amplicon_low, amplicon_high;
    int index;
  };
  
  vector<PrimerPair> pairs;
  unordered_map<string, int> pair_id_to_idx;
  unordered_map<string, vector<int>> locus_to_pairs;
  unordered_map<string, int> locus_to_index;
  vector<string> loci;
  
  int amplicon_len_max=0;
  int amplicon_high_max=0;
  int amplicon_low_min=100000;
  
  for(int i = 0; i < n; i++) {
    string locus = as<string>(str_loci_all[i]);
    string id5 = as<string>(ids_5[i]);
    string id3 = as<string>(ids_3[i]);
    
    PrimerPair pair;
    pair.pair_id = locus + "$" + id5 + "&" + locus + "$" + id3;
    pair.locus = locus;
    pair.id5 = id5;
    pair.id3 = id3;
    pair.weight = weights[i];
    pair.distance = dis_5[i] + dis_3[i];
    pair.amplicon_low = amplicon_size_low[i];
    pair.amplicon_high = amplicon_size_high[i];
    pair.index = pairs.size();
    
    pair.temp5 = (temp_5[i]);
    pair.temp3 = (temp_3[i]);
    
    pair.primers5 = as<string>(primers_5[i]);
    pair.primers3 = as<string>(primers_3[i]);
    
    pairs.push_back(pair);
    pair_id_to_idx[pair.pair_id] = pair.index;
    locus_to_pairs[locus].push_back(pair.index);
    
    if(locus_to_index.find(locus) == locus_to_index.end()) {
      locus_to_index[locus] = loci.size();
      loci.push_back(locus);
    }
    
    if(amplicon_size_high[i] > amplicon_high_max){
      amplicon_high_max=amplicon_size_high[i];
    }
    if(amplicon_size_low[i] < amplicon_low_min){
      amplicon_low_min=amplicon_size_low[i];
    }
    
  }
  
  amplicon_len_max=amplicon_high_max-amplicon_low_min+1;
  
  
  int num_pairs = pairs.size();
  int num_loci = loci.size();
  
  // 2. 建立兼容性集合
  set<pair<int, int>> compatible_set;
  
  for(int i = 0; i < compatible_pairs.size(); i++) {
    string comp_str = as<string>(compatible_pairs[i]);
    size_t pos = comp_str.find("#");
    if(pos != string::npos) {
      string pair1 = comp_str.substr(0, pos);
      string pair2 = comp_str.substr(pos + 1);
      
      if(pair_id_to_idx.find(pair1) != pair_id_to_idx.end() &&
         pair_id_to_idx.find(pair2) != pair_id_to_idx.end()) {
        int idx1 = pair_id_to_idx[pair1];
        int idx2 = pair_id_to_idx[pair2];
        
        if(idx1 < idx2) {
          compatible_set.insert(make_pair(idx1, idx2));
        } else {
          compatible_set.insert(make_pair(idx2, idx1));
        }
      }
    }
  }
  
  // 3. 生成ILP模型
  string lp_filename = "fixed_multiplex_model2_locus_D_correct_binary.lp";
  ofstream lp(lp_filename);
  lp << fixed << setprecision(0);
  
  // 3.1 目标函数
  lp << "Maximize\n";
  lp << "obj: ";
  
  vector<string> Y_vars;  // Y_i (引物对选中变量)
  vector<string> C_vars;  // C_{k,t} (位点k分配到通道t，二进制)
  vector<string> D_vars;  // D_{k,l} 表示位点k和l在同一通道
  vector<string> S_vars;  // S_{i,j,t} 扩增子大小关系（对于通道t）
  unordered_set<string> Z_vars; // Z_k 表示位点k未被选中
  vector<string> E_vars;  // E_{k,l,t} 辅助变量
  
  // 目标函数: 最大化总权重
  for(int i = 0; i < num_pairs; i++) {
    string Y_var = "Y" + to_string(i);
    Y_vars.push_back(Y_var);
    
    double coeff = pairs[i].weight;
    lp << coeff << " " << Y_var;
    
    if(i != num_pairs-1) {
      lp << " + ";
      if((i+1) % 5 == 0) {
        lp << "\n      ";
      }
    }
  }
  lp << "\n\n";
  
  // 约束部分
  lp << "Subject To\n";
  
  // 约束0: 目标函数下界约束（可选）
  if(lower_bound.size() > 0 && lower_bound[0] > 0) {
    for(int i = 0; i < num_pairs; i++) {
      string Y_var = "Y" + to_string(i);
      lp << pairs[i].weight << " " << Y_var;
      
      if(i != num_pairs-1) {
        lp << " + ";
        if((i+1) % 5 == 0) {
          lp << "\n      ";
        }
      }
    }
    lp << " >= " << lower_bound[0] << "\n\n";
  }
  
  //这个约束用来提升求解速度
  //约束0.1，总的扩增子的长度+gap应该小于panel的长度+gap
  for(int i = 0; i < num_pairs; i++) {
    string Y_var = "Y" + to_string(i);
    lp << pairs[i].amplicon_high-pairs[i].amplicon_low+1+ size_gap<< " " << Y_var;
    
    if(i != num_pairs-1) {
      lp << " + ";
      if((i+1) % 5 == 0) {
        lp << "\n      ";
      }
    }
  }
  lp << " <= " << amplicon_len_max*num_channels+ size_gap*num_channels+size_gap<< "\n\n";

  
  
  
  
  
  
  int constr_count = 0;
  
  // 创建C_{k,t}变量（每个位点对应每个通道都有一个二进制变量）
  vector<vector<string>> C_vars_matrix(num_loci, vector<string>(num_channels));
  for(int k = 0; k < num_loci; k++) {
    for(int t = 0; t < num_channels; t++) {
      string C_var = "C" + to_string(k) + "_" + to_string(t);
      C_vars_matrix[k][t] = C_var;
      C_vars.push_back(C_var);
      
      // C_{k,t} 是二进制变量
      lp << C_var << " >= 0\n";
      lp << C_var << " <= 1\n";
      constr_count += 2;
    }
  }
  
  // 约束1: 每个位点最多只能分配到一个通道
  for(int k = 0; k < num_loci; k++) {
    for(int t = 0; t < num_channels; t++) {
      lp << C_vars_matrix[k][t];
      if(t != num_channels-1) lp << " + ";
    }
    lp << " <= 1\n";
    constr_count++;
  }
  
  
  
  
  //    string D_var = "D" + to_string(k) + "_" + to_string(l);
  //    这个约束用来加速
  // 每个位点的引物的集合Y1+Y2+Y3+...>= D1_k
  for(const auto& entry1 : locus_to_pairs) {
    
    int k=locus_to_index.at(entry1.first);
    
    for(const auto& entry2 : locus_to_pairs) {
      int l=locus_to_index.at(entry2.first);
      
      if(k>=l) continue;
      
      string D_var = "D" + to_string(k) + "_" + to_string(l);
      

      
      
      const vector<int>& pair_indices1 = entry1.second;
      
      if(pair_indices1.size() > 1) {
        for(size_t i = 0; i < pair_indices1.size(); i++) {
          int idx = pair_indices1[i];
          lp << "Y" << idx;
          if(i != pair_indices1.size()-1) lp << " + ";
        }
        lp << " - "<< D_var <<" >=0\n";
        constr_count++;
      }
      
      
      const vector<int>& pair_indices2 = entry2.second;
      
      if(pair_indices2.size() > 1) {
        for(size_t i = 0; i < pair_indices2.size(); i++) {
          int idx = pair_indices2[i];
          lp << "Y" << idx;
          if(i != pair_indices2.size()-1) lp << " + ";
        }
        lp << " - "<< D_var <<" >=0\n";
        constr_count++;
      }
      
    
    }
    
  }
  
  
  
  //   
      
      
      
      

  
  
  
  // 约束2: 如果位点k的任意引物对被选中，那么该位点必须被分配到一个通道
  // Y_i ≤ ∑_t C_{k,t} （对于位点k的每个引物对i）
  for(int k = 0; k < num_loci; k++) {
    string locus = loci[k];
    const vector<int>& pair_indices = locus_to_pairs[locus];
    
    for(int idx : pair_indices) {
      string Y_var = "Y" + to_string(idx);
      
      lp << Y_var;
      for(int t = 0; t < num_channels; t++) {
        lp << " - " << C_vars_matrix[k][t];
      }
      lp << " <= 0\n";
      constr_count++;
    }
  }
  
  // 约束3: 如果位点k被分配到通道t，那么该位点必须有一个引物对被选中
  // C_{k,t} ≤ ∑_{i∈I(k)} Y_i
  // 因为每个位点最多只能选一个引物对（约束4），所以∑_{i∈I(k)} Y_i ≤ 1
  for(int k = 0; k < num_loci; k++) {
    string locus = loci[k];
    const vector<int>& pair_indices = locus_to_pairs[locus];
    
    for(int t = 0; t < num_channels; t++) {
      string C_var = C_vars_matrix[k][t];
      
      // C_{k,t} ≤ Y_{i1} + Y_{i2} + ... （所有属于位点k的引物对）
      lp << C_var;
      for(int idx : pair_indices) {
        string Y_var = "Y" + to_string(idx);
        lp << " - " << Y_var;
      }
      lp << " <= 0\n";
      constr_count++;
    }
  }
  
  // 约束4: 同一引物位点最多选一个引物对
  for(const auto& entry : locus_to_pairs) {
    const vector<int>& pair_indices = entry.second;
    
    if(pair_indices.size() > 1) {
      for(size_t i = 0; i < pair_indices.size(); i++) {
        int idx = pair_indices[i];
        lp << "Y" << idx;
        if(i != pair_indices.size()-1) lp << " + ";
      }
      lp << " <= 1\n";
      constr_count++;
    }
  }
  
  // 约束5: 定义每个位点的Z_k变量 (Z_k = 1 - ∑_t C_{k,t})
  // Z_k = 1 表示位点k未被选中
  for(int k = 0; k < num_loci; k++) {
    string Z_var = "Z" + to_string(k);
    Z_vars.insert(Z_var);
    
    lp << Z_var;
    for(int t = 0; t < num_channels; t++) {
      lp << " + " << C_vars_matrix[k][t];
    }
    lp << " = 1\n";
    constr_count++;
  }
  
  // 约束6: 定义位点之间的D_{k,l}变量
  // D_{k,l} = 1 当且仅当 存在某个通道t使得 C_{k,t} = C_{l,t} = 1
  for(int k = 0; k < num_loci; k++) {
    for(int l = k+1; l < num_loci; l++) {
      string D_var = "D" + to_string(k) + "_" + to_string(l);
      D_vars.push_back(D_var);
      
      // 对于每个通道t，定义辅助变量E_{k,l,t} = C_{k,t} * C_{l,t}
      for(int t = 0; t < num_channels; t++) {
        string E_var = "E" + to_string(k) + "_" + to_string(l) + "_" + to_string(t);
        E_vars.push_back(E_var);
        string C_kt = C_vars_matrix[k][t];
        string C_lt = C_vars_matrix[l][t];
        
        // E_{k,l,t} ≤ C_{k,t}
        lp << E_var << " - " << C_kt << " <= 0\n";
        // E_{k,l,t} ≤ C_{l,t}
        lp << E_var << " - " << C_lt << " <= 0\n";
        // E_{k,l,t} ≥ C_{k,t} + C_{l,t} - 1
        lp << E_var << " - " << C_kt << " - " << C_lt << " >= -1\n";
        constr_count += 3;
        
        // D_{k,l} ≥ E_{k,l,t} 对于所有t
        lp << D_var << " - " << E_var << " >= 0\n";
        constr_count++;
      }
      
      // D_{k,l} ≤ ∑_t E_{k,l,t}
      lp << D_var;
      for(int t = 0; t < num_channels; t++) {
        string E_var = "E" + to_string(k) + "_" + to_string(l) + "_" + to_string(t);
        lp << " - " << E_var;
      }
      lp << " <= 0\n";
      constr_count++;
    }
  }
  /*
  // 约束7: 不兼容引物对不能共存
  for(int i = 0; i < num_pairs; i++) {
    for(int j = i+1; j < num_pairs; j++) {
      if(pairs[i].locus == pairs[j].locus) continue;
      
      bool is_compatible = (compatible_set.find(make_pair(i, j)) != compatible_set.end());
      
      if(is_compatible) continue;
      
      string Y_i = "Y" + to_string(i);
      string Y_j = "Y" + to_string(j);
      
      lp << Y_i << " + " << Y_j << " <= 1\n";
      constr_count++;
    }
  }
  */
  
  
  // 10 Y102 + Y134+Y135+Y151+Y152+Y168+Y169+Y186+Y198+Y199 <= 10
  // 10 Y103 + Y134+Y135+Y151+Y152+Y168+Y169+Y185+Y186+Y198+Y199 <= 10
  // 10 Y104 + Y134+Y135+Y151+Y152+Y168+Y169+Y185+Y186+Y198+Y199 <= 10
  // 10 Y105 + Y134+Y135+Y151+Y152+Y168+Y169+Y185+Y186+Y198+Y199 <= 10
  // 10 Y106 + Y134+Y135+Y151+Y152+Y168+Y169+Y185+Y186+Y198+Y200 <= 10
  // 可以合并成 8 Y102 +8 Y103 +8 Y104 + 8 Y105 +  8 Y106 +Y134+Y135+Y151+Y152+Y168+Y169+Y185+Y186+Y198+Y199 <= 8
  // 1 Y103 +Y199 <= 1
  // 2 Y104 +Y185+Y199 <= 2
  // 2 Y105 +Y185+Y199 <= 2
  // 1 Y105 +Y199 <= 1
  // 1 Y106 +Y200 <= 1
  
  // 按位点顺序处理不兼容约束
  // 首先，为每个位点收集所有不兼容关系
  // unordered_map<string, unordered_map<int, set<int>>> locus_incompatible_map; // 位点 -> {y_k -> {不兼容的y_j集合}}
  // std::cout<<"locus: "<<locus<<" get incompatible_sets"<<std::endl;
  // 
  // // 收集所有不兼容关系
  // for(int i = 0; i < num_pairs; i++) {
  //   for(int j = i+1; j < num_pairs; j++) {
  //     if(pairs[i].locus == pairs[j].locus) continue;
  //     
  //     bool is_compatible = (compatible_set.find(make_pair(i, j)) != compatible_set.end());
  //     if(is_compatible) continue;
  //     
  //     // 两个引物对不兼容
  //     string locus_i = pairs[i].locus;
  //     locus_incompatible_map[locus_i][i].insert(j);
  //     
  //     string locus_j = pairs[j].locus;
  //     locus_incompatible_map[locus_j][j].insert(i);
  //   }
  // }
  
  // 按位点顺序处理不兼容约束
  unordered_set<string> covered_constraints; // 记录已经覆盖的不兼容关系
  
  for(const auto& locus_entry : locus_to_pairs) {
    const string& locus = locus_entry.first;
    const vector<int>& pair_indices = locus_entry.second;
    std::cout<<"locus: "<<locus<<std::endl;
    
    // 如果这个位点只有1个引物对
    if(pair_indices.size() <= 1) {
      for(int y_k : pair_indices) {
        set<int> incompatible_set;
        
        // 收集所有不兼容的引物对，排除已覆盖的
        for(int j = 0; j < num_pairs; j++) {
          if(pairs[y_k].locus == pairs[j].locus) continue;
          
          bool is_compatible = (compatible_set.find(make_pair(min(y_k, j), max(y_k, j))) != compatible_set.end());
          if(is_compatible) continue;
          
          // 检查这个约束是否已经被覆盖
          string constraint_key = to_string(min(y_k, j)) + "_" + to_string(max(y_k, j));
          if(covered_constraints.find(constraint_key) == covered_constraints.end()) {
            incompatible_set.insert(j);
          }
        }
        
        if(!incompatible_set.empty()) {
          int n = incompatible_set.size();
          string concat_str = "";
          for(int y_j : incompatible_set) {
            concat_str += "Y" + to_string(y_j) + " + ";
            
            // 标记这个约束为已覆盖
            string key1 = to_string(y_k) + "_" + to_string(y_j);
            string key2 = to_string(y_j) + "_" + to_string(y_k);
            covered_constraints.insert(key1);
            covered_constraints.insert(key2);
          }
          if(!concat_str.empty()) {
            concat_str.pop_back(); // 移除最后一个" +"
            concat_str.pop_back();
          }
          
          lp << n << " Y" << y_k << " + " << concat_str << " <= " << n << "\n";
          constr_count++;
        }
      }
      continue;
    }
    std::cout<<"locus: "<<locus<<" get incompatible_sets"<<std::endl;
    
    // 位点有多个引物对，收集每个引物对的不兼容关系
    vector<set<int>> incompatible_sets(pair_indices.size());
    
    for(size_t idx = 0; idx < pair_indices.size(); idx++) {
      int y_k = pair_indices[idx];
      
      // 收集y_k的所有不兼容关系，排除已覆盖的
      for(int j = 0; j < num_pairs; j++) {
        if(pairs[y_k].locus == pairs[j].locus) continue;
        
        bool is_compatible = (compatible_set.find(make_pair(min(y_k, j), max(y_k, j))) != compatible_set.end());
        if(is_compatible) continue;
        
        // 检查这个约束是否已经被覆盖
        string constraint_key = to_string(min(y_k, j)) + "_" + to_string(max(y_k, j));
        if(covered_constraints.find(constraint_key) == covered_constraints.end()) {
          incompatible_sets[idx].insert(j);
        }
      }
    }
    
    
    //std::cout<<"locus: "<<locus<<" get simimarity matrix"<<std::endl;
    
    
    // 优先两两聚类, 后续的层次聚类要求高
    // 在这个位点内部进行聚类
    // 在这个位点内部进行聚类
    int n = pair_indices.size();
    
    // 计算相似度矩阵
    std::cout << "locus: " << locus << " hierarchical clustering" << std::endl;
    
    const int MIN_SHARED_ELEMENTS = max(1, num_pairs / 100);
    const double BASE_COVERAGE_RATIO = 0.3;
    const double MAX_ADDITIONAL = 0.9;
    const double SIZE_PENALTY_RATE = 0.2;
    
    // 新的数据结构定义
    struct QueueItem {
      double coverage;
      int intersection;
      int i, j;
      
      // 自定义比较函数：coverage大的在前（降序）
      bool operator<(const QueueItem& other) const {
        if (std::fabs(coverage - other.coverage) > 1e-9)
          return coverage > other.coverage;  // 降序
        if (i != other.i) return i < other.i;
        return j < other.j;
      }
    };
    
    struct PairHash {
      std::size_t operator()(const std::pair<int, int>& p) const noexcept {
        return std::hash<int>()(p.first) ^ (std::hash<int>()(p.second) << 1);
      }
    };
    
    set<QueueItem> cluster_queue;  // 自动按coverage降序排序
    unordered_map<pair<int, int>, set<QueueItem>::iterator, PairHash> queue_index;
    
    // 辅助函数：更新或插入队列项
    auto updateQueue = [&](int i, int j, double coverage, int intersection) {
      // 统一键值：确保(i,j)和(j,i)被视为同一个pair
      int id1 = min(i, j);
      int id2 = max(i, j);
      auto key = make_pair(id1, id2);
      
      // 如果已经存在，先删除旧的
      if (queue_index.find(key) != queue_index.end()) {
        cluster_queue.erase(queue_index[key]);
        queue_index.erase(key);
      }
      
      // 只有可能合并的才加入队列
      if (intersection >= MIN_SHARED_ELEMENTS) {
        QueueItem new_item{coverage, intersection, i, j};
        auto result = cluster_queue.insert(new_item);
        if (result.second) {  // 插入成功
          queue_index[key] = result.first;
        }
      }
    };
    
    auto cleanupClusterFromQueue = [&](int cluster_id) {
      auto it = queue_index.begin();
      while (it != queue_index.end()) {
        if (it->first.first == cluster_id || it->first.second == cluster_id) {
          cluster_queue.erase(it->second);
          it = queue_index.erase(it);
        } else {
          ++it;
        }
      }
    };
    
    
    // 1. 初始化：每个元素一个cluster
    vector<vector<int>> clusters(n);
    vector<set<int>> cluster_centers(n);
    vector<int> cluster_sizes(n, 1);
    
    for(int i = 0; i < n; i++) {
      clusters[i] = {i};
      cluster_centers[i] = incompatible_sets[i];
    }
    
    // 2. 预计算相似度并初始化队列
    for(int i = 0; i < n; i++) {
      for(int j = i+1; j < n; j++) {
        const set<int>& set1 = incompatible_sets[i];
        const set<int>& set2 = incompatible_sets[j];
        
        if(set1.empty() || set2.empty()) continue;
        
        // 计算交集大小
        int intersection = 0;
        auto it1 = set1.begin();
        auto it2 = set2.begin();
        while(it1 != set1.end() && it2 != set2.end()) {
          if(*it1 < *it2) ++it1;
          else if(*it2 < *it1) ++it2;
          else {
            intersection++;
            ++it1;
            ++it2;
          }
        }
        
        if(intersection < MIN_SHARED_ELEMENTS) continue;
        
        // 选择使用min还是max计算覆盖率（根据需求）
        // 选项A：使用min（更严格）
        // int min_size = min(set1.size(), set2.size());
        // double coverage = (double)intersection / min_size;
        
        // 选项B：使用max（更宽松）
        // 使用较大的cluster对应的center的size
        int min_size_new = (set1.size() <= set2.size()) ? set1.size() : set2.size();
        double coverage = (double)intersection / min_size_new;
        
        if(coverage >= BASE_COVERAGE_RATIO) {
          updateQueue(i, j, coverage, intersection);
        }
        
      }
    }
    
    // 3. 维护cluster有效性标志
    vector<bool> cluster_active(n, true);
    
    // 4. 层次聚类主循环
    while(!cluster_queue.empty()) {
      // 获取最佳pair
      auto best_it = cluster_queue.begin();
      QueueItem best_item = *best_it;
      
      double coverage = best_item.coverage;
      int intersection = best_item.intersection;
      int i = best_item.i;
      int j = best_item.j;
      
      // 从队列中移除
      auto key = make_pair(min(i, j), max(i, j));
      cluster_queue.erase(best_it);
      queue_index.erase(key);
      
      // 检查cluster是否仍然有效
      if(!cluster_active[i] || !cluster_active[j]) continue;
      
      // 计算动态阈值（保持你的计算方式）
      double dynamic_threshold = (cluster_sizes[i] > 1 || cluster_sizes[j] > 1) ? 
      0.9 : BASE_COVERAGE_RATIO;
      
      if(coverage < dynamic_threshold) continue;
      
      // 检查合并后的交集
      set<int> new_center;
      set_intersection(cluster_centers[i].begin(), cluster_centers[i].end(),
                       cluster_centers[j].begin(), cluster_centers[j].end(),
                       inserter(new_center, new_center.begin()));
      
      if(new_center.size() < MIN_SHARED_ELEMENTS) continue;
      
      // 合并cluster j 到 i
      clusters[i].insert(clusters[i].end(), clusters[j].begin(), clusters[j].end());
      cluster_centers[i] = new_center;
      cluster_sizes[i] += cluster_sizes[j];
      cluster_active[j] = false;
      
      cleanupClusterFromQueue(j);  // 清理所有包含j的记录
      
      // 可选：清理i的旧记录（因为i的center已更新）
      cleanupClusterFromQueue(i);  // 先清理所有i的旧记录
      
      
      // 重新计算i与其他clusters的相似度
      for(int k = 0; k < n; k++) {
        if(k == i || !cluster_active[k]) continue;
        
        const set<int>& center_k = cluster_centers[k];
        if(center_k.empty() || new_center.empty()) continue;
        
        // 计算新的交集
        int new_intersection = 0;
        auto it_i = new_center.begin();
        auto it_k = center_k.begin();
        while(it_i != new_center.end() && it_k != center_k.end()) {
          if(*it_i < *it_k) ++it_i;
          else if(*it_k < *it_i) ++it_k;
          else {
            new_intersection++;
            ++it_i;
            ++it_k;
          }
        }
        
        // 计算新的覆盖度（保持你的计算方式）
        if (new_intersection >= MIN_SHARED_ELEMENTS) {
          // 使用较大的cluster对应的center的size
          int max_size_new = (cluster_sizes[i] >= cluster_sizes[k]) ? 
          new_center.size() : center_k.size();
          double new_coverage = (double)new_intersection / max_size_new;
          
          // 动态阈值（保持你的计算方式）
          double dynamic_threshold_k = (cluster_sizes[i] > 1 || cluster_sizes[k] > 1) ? 
          0.9 : BASE_COVERAGE_RATIO;
          
          // 总是更新队列，无论是否满足阈值
          updateQueue(i, k, new_coverage, new_intersection);
        } else {
          // 交集太小，确保从队列中删除（如果存在）
          updateQueue(i, k, 0.0, new_intersection);
        }
      }
      
      // 也需要更新其他cluster与已删除cluster j的相关记录
      // 这些记录在后续获取时会因为cluster_active[j]=false而被跳过
    }
    
    // 5. 收集结果
    vector<vector<int>> valid_clusters;
    vector<set<int>> valid_centers;
    
    for(int i = 0; i < n; i++) {
      if(cluster_active[i] && !clusters[i].empty()) {
        valid_clusters.push_back(clusters[i]);
        valid_centers.push_back(cluster_centers[i]);
      }
    }
    
    clusters = valid_clusters;
    cluster_centers = valid_centers;
    
    
    
    
    
    
    std::cout<<"locus: "<<locus<<" get constraints"<<std::endl;
    
    // 为每个cluster生成约束
    for(const auto& cluster : clusters) {
      if(cluster.empty()) continue;
      
      // 找出共享的不兼容集合
      set<int> shared_incompatible;
      if(!cluster.empty()) {
        int first_idx = cluster[0];
        shared_incompatible = incompatible_sets[first_idx];
        
        for(size_t i = 1; i < cluster.size(); i++) {
          int idx = cluster[i];
          const set<int>& current_set = incompatible_sets[idx];
          
          set<int> intersection;
          set_intersection(shared_incompatible.begin(), shared_incompatible.end(),
                           current_set.begin(), current_set.end(),
                           inserter(intersection, intersection.begin()));
          shared_incompatible = intersection;
        }
      }
      
      // 计算特有集合
      vector<set<int>> specific_sets(cluster.size());
      for(size_t i = 0; i < cluster.size(); i++) {
        int idx = cluster[i];
        const set<int>& all_incompatible = incompatible_sets[idx];
        
        set<int> specific;
        set_difference(all_incompatible.begin(), all_incompatible.end(),
                       shared_incompatible.begin(), shared_incompatible.end(),
                       inserter(specific, specific.begin()));
        specific_sets[i] = specific;
      }
      
      // 生成共享约束
      if(cluster.size() > 1 && !shared_incompatible.empty()) {
        int n_shared = shared_incompatible.size();
        
        bool first = true;
        lp << n_shared;
        for(int idx : cluster) {
          int y_k = pair_indices[idx];
          if(first) {
            lp << " Y" << y_k;
            first = false;
          } else {
            lp << " + " << n_shared << " Y" << y_k;
          }
        }
        
        for(int y_j : shared_incompatible) {
          lp << " + Y" << y_j;
        }
        
        lp << " <= " << n_shared << "\n";
        constr_count++;
        
        // 标记所有共享的不兼容关系为已覆盖
        for(int idx : cluster) {
          int y_k = pair_indices[idx];
          for(int y_j : shared_incompatible) {
            string key1 = to_string(y_k) + "_" + to_string(y_j);
            string key2 = to_string(y_j) + "_" + to_string(y_k);
            covered_constraints.insert(key1);
            covered_constraints.insert(key2);
          }
        }
      }
      
      // 生成特有约束
      for(size_t i = 0; i < cluster.size(); i++) {
        int idx = cluster[i];
        int y_k = pair_indices[idx];
        const set<int>& specific_set = specific_sets[i];
        
        if(!specific_set.empty()) {
          int n_specific = specific_set.size();
          
          lp << n_specific << " Y" << y_k;
          for(int y_j : specific_set) {
            lp << " + Y" << y_j;
          }
          lp << " <= " << n_specific << "\n";
          constr_count++;
          
          // 标记特有约束为已覆盖
          for(int y_j : specific_set) {
            string key1 = to_string(y_k) + "_" + to_string(y_j);
            string key2 = to_string(y_j) + "_" + to_string(y_k);
            covered_constraints.insert(key1);
            covered_constraints.insert(key2);
          }
        }
      }
      
      // 处理单个Y的情况
      if(cluster.size() == 1) {
        int idx = cluster[0];
        int y_k = pair_indices[idx];
        const set<int>& all_incompatible = incompatible_sets[idx];
        
        if(!all_incompatible.empty()) {
          int n_all = all_incompatible.size();
          lp << n_all << " Y" << y_k;
          for(int y_j : all_incompatible) {
            lp << " + Y" << y_j;
          }
          lp << " <= " << n_all << "\n";
          constr_count++;
          
          // 标记为已覆盖
          for(int y_j : all_incompatible) {
            string key1 = to_string(y_k) + "_" + to_string(y_j);
            string key2 = to_string(y_j) + "_" + to_string(y_k);
            covered_constraints.insert(key1);
            covered_constraints.insert(key2);
          }
        }
      }
    }
  }
  
  
  
  
  
  
  
  // Y54 + Y102 + Y104 + Y108 + Y109 + Y110  + Y116 + D3_6 <= 2
  // Y55 + Y102 + Y104 + Y108 + Y109 + Y110 + Y113  + D3_6 <= 2
  // 可以合并成
  // Y54 + Y55 + Y102 + Y104 + Y108 + Y109 + Y110 + D3_6 <= 2
  // 
  // Y54 + Y116 + D3_6 <= 2
  // Y55 + Y113 + D3_6 <= 2
  // 约束8-9: 同一通道内兼容引物对的扩增子重叠约束
  // 使用D_{k,l}变量（位点k和l是否在同一通道）
  // 约束8: 兼容但可能重叠的引物对不能在同一通道
  // 约束8: 兼容但可能重叠的引物对不能在同一通道
  // 使用映射来合并约束
  // 约束8-9: 同一通道内兼容引物对的扩增子重叠约束
  // 优化版本：基于不兼容模式相似度的聚类
  // 约束8-9: 同一通道内兼容引物对的扩增子重叠约束
  // 优化版本：基于相似度的聚类合并
  unordered_map<string, unordered_map<int, set<int>>> locus_pair_to_y_pairs;
  
  std::cout<<"D_var: "<<" map to D_var"<<std::endl;
  
  // 1. 收集所有相关Y对
  for(int i = 0; i < num_pairs; i++) {
    for(int j = i+1; j < num_pairs; j++) {
      if(pairs[i].locus == pairs[j].locus) continue;
      
      bool is_compatible = (compatible_set.find(make_pair(i, j)) != compatible_set.end());
      if(!is_compatible) continue;
      
      bool could_overlap = !(pairs[i].amplicon_high + size_gap <= pairs[j].amplicon_low ||
                             pairs[j].amplicon_high + size_gap <= pairs[i].amplicon_low);
      if(!could_overlap) continue;
      
      int locus_idx_i = locus_to_index[pairs[i].locus];
      int locus_idx_j = locus_to_index[pairs[j].locus];
      int k = min(locus_idx_i, locus_idx_j);
      int l = max(locus_idx_i, locus_idx_j);
      
      string D_var = "D" + to_string(k) + "_" + to_string(l);
      
      int y_k, y_l;
      if(locus_to_index[pairs[i].locus] == k) {
        y_k = i;
        y_l = j;
      } else {
        y_k = j;
        y_l = i;
      }
      
      locus_pair_to_y_pairs[D_var][y_k].insert(y_l);
    }
  }
  
  // 2. 为每个D_var生成优化约束
  for(const auto& d_entry : locus_pair_to_y_pairs) {
    const string& D_var = d_entry.first;
    const auto& y_map = d_entry.second;
    
    if(y_map.empty()) continue;
    
    // 解析位点索引
    size_t underscore1 = D_var.find('_');
    size_t underscore2 = D_var.find('_', underscore1 + 1);
    int k = stoi(D_var.substr(1, underscore1 - 1));
    
    // 找出属于位点k的所有Y_k
    vector<int> yk_list;
    for(const auto& entry : y_map) {
      int y_k = entry.first;
      if(locus_to_index[pairs[y_k].locus] == k) {
        yk_list.push_back(y_k);
      }
    }
    
    if(yk_list.empty()) continue;
    
    // 3. 对Y_k进行聚类（基于不兼容模式的相似度）
    // 每个cluster包含一组相似的Y_k和它们共享的不兼容Y_l
    
    // 3.1 构建相似度矩阵
    int n = yk_list.size();
    vector<vector<double>> similarity_matrix(n, vector<double>(n, 0.0));
    vector<set<int>> incompatible_sets(n);
    
    for(int i = 0; i < n; i++) {
      incompatible_sets[i] = y_map.at(yk_list[i]);
    }
    
    std::cout<<"D_var: "<<D_var<<" get clusters"<<std::endl;
    
    // 计算相似度（Jaccard相似系数）
    for(int i = 0; i < n; i++) {
      for(int j = i+1; j < n; j++) {
        const set<int>& set1 = incompatible_sets[i];
        const set<int>& set2 = incompatible_sets[j];
        
        // 计算交集和并集
        set<int> intersection, union_set;
        set_intersection(set1.begin(), set1.end(),
                         set2.begin(), set2.end(),
                         inserter(intersection, intersection.begin()));
        set_union(set1.begin(), set1.end(),
                  set2.begin(), set2.end(),
                  inserter(union_set, union_set.begin()));
        
        // Jaccard相似度 = |交集| / |并集|
        if(!union_set.empty()) {
          similarity_matrix[i][j] = similarity_matrix[j][i] = 
            (double)intersection.size() / union_set.size();
        }
      }
    }
    
    // 3.2 约束8-9的官方层次聚类算法（使用priority queue）
    std::cout << "D_var: " << D_var << " hierarchical clustering (priority queue)" << std::endl;
    
    const int MIN_SHARED_ELEMENTS = 1;
    const double BASE_COVERAGE_RATIO = 0.3;
    const double MAX_ADDITIONAL = 0.9;
    const double SIZE_PENALTY_RATE = 0.2;
    
    // 1. 初始化：每个Y_k作为一个单独的cluster
    vector<vector<int>> clusters(n);
    vector<set<int>> cluster_centers(n);
    vector<int> cluster_sizes(n, 1);
    
    for(int i = 0; i < n; i++) {
      clusters[i] = {i};
      cluster_centers[i] = incompatible_sets[i];
    }
    
    // 新的数据结构定义
    struct QueueItem {
      double coverage;
      int intersection;
      int i, j;
      
      // 自定义比较函数：coverage大的在前（降序）
      bool operator<(const QueueItem& other) const {
        if (std::fabs(coverage - other.coverage) > 1e-9)
          return coverage > other.coverage;  // 降序
        if (i != other.i) return i < other.i;
        return j < other.j;
      }
    };
    
    struct PairHash {
      std::size_t operator()(const std::pair<int, int>& p) const noexcept {
        return std::hash<int>()(p.first) ^ (std::hash<int>()(p.second) << 1);
      }
    };
    
    set<QueueItem> cluster_queue;  // 自动按coverage降序排序
    unordered_map<pair<int, int>, set<QueueItem>::iterator, PairHash> queue_index;
    
    // 辅助函数：更新或插入队列项
    auto updateQueue = [&](int i, int j, double coverage, int intersection) {
      // 统一键值：确保(i,j)和(j,i)被视为同一个pair
      int id1 = min(i, j);
      int id2 = max(i, j);
      auto key = make_pair(id1, id2);
      
      // 如果已经存在，先删除旧的
      if (queue_index.find(key) != queue_index.end()) {
        cluster_queue.erase(queue_index[key]);
        queue_index.erase(key);
      }
      
      // 只有可能合并的才加入队列
      if (intersection >= MIN_SHARED_ELEMENTS) {
        QueueItem new_item{coverage, intersection, i, j};
        auto result = cluster_queue.insert(new_item);
        if (result.second) {  // 插入成功
          queue_index[key] = result.first;
        }
      }
    };
    
    // 清理函数
    auto cleanupClusterFromQueue = [&](int cluster_id) {
      auto it = queue_index.begin();
      while (it != queue_index.end()) {
        if (it->first.first == cluster_id || it->first.second == cluster_id) {
          cluster_queue.erase(it->second);
          it = queue_index.erase(it);
        } else {
          ++it;
        }
      }
    };
    
    // 3. 预计算相似度并初始化队列
    for(int i = 0; i < n; i++) {
      for(int j = i+1; j < n; j++) {
        const set<int>& set1 = incompatible_sets[i];
        const set<int>& set2 = incompatible_sets[j];
        
        if(set1.empty() || set2.empty()) continue;
        
        // 计算交集大小
        int intersection = 0;
        auto it1 = set1.begin();
        auto it2 = set2.begin();
        while(it1 != set1.end() && it2 != set2.end()) {
          if(*it1 < *it2) ++it1;
          else if(*it2 < *it1) ++it2;
          else {
            intersection++;
            ++it1;
            ++it2;
          }
        }
        
        if(intersection < MIN_SHARED_ELEMENTS) continue;
        
        // 计算覆盖率：使用max
        int min_size = min(set1.size(), set2.size());
        double coverage = (double)intersection / min_size;
        
        if(coverage >= BASE_COVERAGE_RATIO) {
          updateQueue(i, j, coverage, intersection);
        }
      }
    }
    
    // 4. 维护cluster有效性标志
    vector<bool> cluster_active(n, true);
    
    // 5. 层次聚类主循环
    while(!cluster_queue.empty()) {
      // 获取最佳pair
      auto best_it = cluster_queue.begin();
      QueueItem best_item = *best_it;
      
      double coverage = best_item.coverage;
      int intersection = best_item.intersection;
      int i = best_item.i;
      int j = best_item.j;
      
      // 从队列中移除
      auto key = make_pair(min(i, j), max(i, j));
      cluster_queue.erase(best_it);
      queue_index.erase(key);
      
      // 检查cluster是否仍然有效
      if(!cluster_active[i] || !cluster_active[j]) continue;
      
      // 动态阈值（使用原来的计算方式）
      //int max_cluster_size = max(cluster_sizes[i], cluster_sizes[j]);
      //double dynamic_threshold = BASE_COVERAGE_RATIO + 
      //  min(MAX_ADDITIONAL, SIZE_PENALTY_RATE * (max_cluster_size - 1));
      


      // 动态阈值（保持你的计算方式）
      double dynamic_threshold = (cluster_sizes[i] > 1 || cluster_sizes[j] > 1) ? 
      0.9 : BASE_COVERAGE_RATIO;
      
      
      
      if(coverage < dynamic_threshold) continue;
      
      // 检查合并后的交集
      set<int> new_center;
      set_intersection(cluster_centers[i].begin(), cluster_centers[i].end(),
                       cluster_centers[j].begin(), cluster_centers[j].end(),
                       inserter(new_center, new_center.begin()));
      
      if(new_center.size() < MIN_SHARED_ELEMENTS) continue;
      
      // 合并cluster j 到 i
      clusters[i].insert(clusters[i].end(), clusters[j].begin(), clusters[j].end());
      cluster_centers[i] = new_center;
      cluster_sizes[i] += cluster_sizes[j];
      cluster_active[j] = false;
      
      // 清理队列中所有包含j的记录
      cleanupClusterFromQueue(j);
      
      // 清理i的旧记录（因为i的center已更新）
      cleanupClusterFromQueue(i);
      
      // 重新计算i与其他clusters的相似度
      for(int k = 0; k < n; k++) {
        if(k == i || !cluster_active[k]) continue;
        
        const set<int>& center_k = cluster_centers[k];
        if(center_k.empty() || new_center.empty()) continue;
        
        // 计算新的交集
        int new_intersection = 0;
        auto it_i = new_center.begin();
        auto it_k = center_k.begin();
        while(it_i != new_center.end() && it_k != center_k.end()) {
          if(*it_i < *it_k) ++it_i;
          else if(*it_k < *it_i) ++it_k;
          else {
            new_intersection++;
            ++it_i;
            ++it_k;
          }
        }
        
        // 计算新的覆盖度
        if (new_intersection >= MIN_SHARED_ELEMENTS) {
          // 使用max计算覆盖率
          //int max_size_new = max(new_center.size(), center_k.size());
          //double new_coverage = (double)new_intersection / max_size_new;
          
          
          
          int max_size_new = (cluster_sizes[i] >= cluster_sizes[k]) ? 
          new_center.size() : center_k.size();
          double new_coverage = (double)new_intersection / max_size_new;
          
          // 动态阈值（保持你的计算方式）
          //double dynamic_threshold_k = (cluster_sizes[i] > 1 || cluster_sizes[k] > 1) ? 
          //0.9 : BASE_COVERAGE_RATIO;
          
          // 总是更新队列，无论是否满足阈值
          updateQueue(i, k, new_coverage, new_intersection);
        } else {
          // 交集太小，确保从队列中删除（如果存在）
          updateQueue(i, k, 0.0, new_intersection);
        }
      }
    }
    
    // 6. 收集有效clusters
    vector<vector<int>> valid_clusters;
    vector<set<int>> valid_centers;
    
    for(int i = 0; i < n; i++) {
      if(cluster_active[i] && !clusters[i].empty()) {
        valid_clusters.push_back(clusters[i]);
        valid_centers.push_back(cluster_centers[i]);
      }
    }
    
    clusters = valid_clusters;
    cluster_centers = valid_centers;
    
    // 确保所有Y_k都被分配
    // for(int i = 0; i < n; i++) {
    //   if(!assigned[i]) {
    //     clusters.push_back({i});
    //   }
    // }
    std::cout<<"D_var: "<<D_var<<" get constraints"<<std::endl;
    
    
    // 4. 为每个cluster生成约束
    for(const auto& cluster : clusters) {
      if(cluster.empty()) continue;
      
      // 4.1 找出cluster中所有Y_k共享的不兼容Y_l
      set<int> shared_y_l;
      if(!cluster.empty()) {
        int first_idx = cluster[0];
        shared_y_l = incompatible_sets[first_idx];
        
        for(size_t i = 1; i < cluster.size(); i++) {
          int idx = cluster[i];
          const set<int>& current_set = incompatible_sets[idx];
          
          set<int> intersection;
          set_intersection(shared_y_l.begin(), shared_y_l.end(),
                           current_set.begin(), current_set.end(),
                           inserter(intersection, intersection.begin()));
          shared_y_l = intersection;
        }
      }
      
      // 4.2 生成共享约束（如果cluster有多个Y_k且有共享的Y_l）
      if(cluster.size() > 1 && !shared_y_l.empty()) {
        // ΣY_k + Σ共享Y_l + D ≤ 2
        bool first = true;
        lp << " ";
        for(int idx : cluster) {
          int y_k = yk_list[idx];
          if(first) {
            lp << "Y" << y_k;
            first = false;
          } else {
            lp << " + Y" << y_k;
          }
        }
        
        for(int y_l : shared_y_l) {
          lp << " + Y" << y_l;
        }
        
        lp << " + " << D_var << " <= 2\n";
        constr_count++;
      }
      
      // 4.3 为cluster中的每个Y_k生成约束
      for(int idx : cluster) {
        int y_k = yk_list[idx];
        const set<int>& all_incompatible = incompatible_sets[idx];
        
        // 确定需要生成约束的Y_l
        set<int> y_l_to_constrain;
        
        if(cluster.size() > 1 && !shared_y_l.empty()) {
          // 对于多个Y_k的cluster，生成特有部分
          set_difference(all_incompatible.begin(), all_incompatible.end(),
                         shared_y_l.begin(), shared_y_l.end(),
                         inserter(y_l_to_constrain, y_l_to_constrain.begin()));
        } else {
          // 单个Y_k或没有共享Y_l，使用所有不兼容的Y_l
          y_l_to_constrain = all_incompatible;
        }
        
        // 将Y_l_to_constrain分组，避免过长的约束
        // 如果Y_l太多，可以进一步分组
        const int MAX_TERMS_PER_CONSTRAINT = 3000; // 可调整
        
        vector<int> y_l_vector(y_l_to_constrain.begin(), y_l_to_constrain.end());
        for(size_t start = 0; start < y_l_vector.size(); start += MAX_TERMS_PER_CONSTRAINT) {
          size_t end = min(start + MAX_TERMS_PER_CONSTRAINT, y_l_vector.size());
          
          lp << " Y" << y_k;
          for(size_t i = start; i < end; i++) {
            lp << " + Y" << y_l_vector[i];
          }
          lp << " + " << D_var << " <= 2\n";
          constr_count++;
        }
      }
    }
    
    // 5. 处理从另一个方向的不兼容关系（如果需要）
    // 上面的代码只处理了从位点k到位点l的关系
    // 如果还需要处理从位点l到位点k的关系，需要类似的代码
  }
  
  
  // 变量类型定义
  lp << "\nBinary\n";
  
  // Y变量
  for(int i = 0; i < num_pairs; i++) {
    lp << "Y" << i << "\n";
  }
  
  // C变量（位点-通道分配变量）
  for(const string& c_var : C_vars) {
    lp << c_var << "\n";
  }
  
  // D变量（位点之间的）
  for(const string& d_var : D_vars) {
    lp << d_var << "\n";
  }
  
  // S变量（引物对之间的）
  for(const string& s_var : S_vars) {
    lp << s_var << "\n";
  }
  
  // Z变量（位点的）
  for(const string& z_var : Z_vars) {
    lp << z_var << "\n";
  }
  
  // E变量（辅助变量）
  for(const string& e_var : E_vars) {
    lp << e_var << "\n";
  }
  
  lp << "\nEnd\n";
  lp.close();
  
  // 统计信息
  int total_incompatible_pairs = 0;
  int total_overlap_pairs = 0;
  
  for(int i = 0; i < num_pairs; i++) {
    for(int j = i+1; j < num_pairs; j++) {
      if(pairs[i].locus == pairs[j].locus) continue;
      
      bool is_compatible = (compatible_set.find(make_pair(i, j)) != compatible_set.end());
      
      if(!is_compatible) {
        total_incompatible_pairs++;
      } else {
        bool could_overlap = !(pairs[i].amplicon_high + size_gap < pairs[j].amplicon_low ||
                               pairs[j].amplicon_high + size_gap < pairs[i].amplicon_low);
        if(could_overlap) {
          total_overlap_pairs++;
        }
      }
    }
  }
  
  // 保存引物信息
  string primers_filename = "primers_info_locus_D_correct_binary.tsv";
  ofstream primers_file(primers_filename);
  primers_file << "pair_index\tpair_id\tlocus\tlocus_index\tid5\tid3\tprimer5\tprimer3\ttemp5\ttemp3\tamplicon_low\tamplicon_high\tweight" << endl;
  
  for(int i = 0; i < num_pairs; i++) {
    primers_file << i << "\t"
                 << pairs[i].pair_id << "\t"
                 << pairs[i].locus << "\t"
                 << locus_to_index[pairs[i].locus] << "\t"
                 << pairs[i].id5 << "\t"
                 << pairs[i].id3 << "\t"
                 << pairs[i].primers5 << "\t"
                 << pairs[i].primers3 << "\t"
                 << pairs[i].temp5 << "\t"
                 << pairs[i].temp3 << "\t"
                 << pairs[i].amplicon_low << "\t"
                 << pairs[i].amplicon_high << "\t"
                 << pairs[i].weight << endl;
  }
  primers_file.close();
  
  // 预览LP文件前1000个字符
  ifstream preview_file(lp_filename);
  string lp_preview(1000, ' ');
  preview_file.read(&lp_preview[0], 1000);
  preview_file.close();
  
  // 返回结果
  List result = List::create(
    Named("primers_info_file") = primers_filename,
    Named("lp_file") = lp_filename,
    Named("num_variables") = Y_vars.size() + C_vars.size() + D_vars.size() + S_vars.size() + 
      Z_vars.size() + E_vars.size(),
      Named("num_constraints") = constr_count,
      Named("num_pairs") = num_pairs,
      Named("num_loci") = num_loci,
      Named("num_channels") = num_channels,
      Named("compatible_pairs") = (int)compatible_set.size(),
      Named("incompatible_pairs") = total_incompatible_pairs,
      Named("overlap_pairs") = total_overlap_pairs,
      Named("model_summary") = List::create(
        Named("Y_vars") = (int)Y_vars.size(),
        Named("C_vars") = (int)C_vars.size(),  // num_loci * num_channels
        Named("D_vars") = (int)D_vars.size(),
        Named("S_vars") = (int)S_vars.size(),
        Named("Z_vars") = (int)Z_vars.size(),
        Named("E_vars") = (int)E_vars.size()
      ),
      Named("lp_preview") = lp_preview + "..."
  );
  
  return result;
}