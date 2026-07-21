#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <set>
#include <map>
#include <cmath>
#include <array>
#include <iostream>
#include <iomanip>
#include <unordered_set>
#include <functional>
#include <queue>
#include <climits> 
#include <stdexcept>
#include <getopt.h>

struct SDFAtom {
    int id;
    double x, y, z;
    std::string element;
    std::string name;
    std::string type;
    std::vector<int> bonds;
    bool inRing;
    bool inAromaticRing;
    int massDifference;
    int charge;
    int stereoParity;
    int hydrogenCount;
    int stereoCount;
    int valence;
    std::string extraInfo;  // 用于存储其他可能的信息   
    SDFAtom() : x(0), y(0), z(0), inRing(false), inAromaticRing(false), massDifference(0),
                charge(0), stereoParity(0), hydrogenCount(0), stereoCount(0), valence(0) {}
};

struct SDFBond {
    int atom1, atom2;
    int type;
    bool rotatable;
    bool inRing;
    bool isAromatic;
    std::string exinfo;  // 用于存储其他可能的信息
    SDFBond() : atom1(-1), atom2(-1), type(0), rotatable(false), inRing(false), isAromatic(false) {}  // 默认构造函数

};

struct SDFMolecule {
    std::string name;
    std::vector<SDFAtom> atoms;
    std::vector<SDFBond> bonds;
    std::vector<std::vector<int>> rings;
    std::string MInfo;  // 存储 M 开头的行信息
    std::string MProp;  // 存储分子属性信息
    bool ringsComputed;  // 标记环是否已计算
};

// 前向声明
void findAllRings(SDFMolecule& mol);

// 安全的字符串到整数转换
int safeStoi(const std::string& str, int defaultValue = 0) {
    try {
        return std::stoi(str);
    } catch (const std::exception& e) {
        std::cerr << "Warning: Failed to convert '" << str << "' to int. Using default value " << defaultValue << std::endl;
        return defaultValue;
    }
}

// 安全的字符串到双精度浮点数转换
double safeStod(const std::string& str, double defaultValue = 0.0) {
    try {
        return std::stod(str);
    } catch (const std::exception& e) {
        std::cerr << "Warning: Failed to convert '" << str << "' to double. Using default value " << defaultValue << std::endl;
        return defaultValue;
    }
}

SDFMolecule parseSDFImpl(std::istream& input) {
    SDFMolecule mol;
    std::string line;
    std::map<std::string, int> elementCount;

    // Read molecule name (first line)
    std::getline(input, mol.name);
    mol.name = mol.name.substr(0, mol.name.find_first_of("\t\r\n"));  // Remove any trailing whitespace

    // Skip next two lines
    std::getline(input, line);
    std::getline(input, line);

    // Read counts line
    std::getline(input, line);
    int atomCount = safeStoi(line.substr(0, 3));
    int bondCount = safeStoi(line.substr(3, 3));

    // Read atoms
    for (int i = 0; i < atomCount; ++i) {
        std::getline(input, line);
        SDFAtom atom;
        atom.id = i;
        atom.x = safeStod(line.substr(0, 10));
        atom.y = safeStod(line.substr(10, 10));
        atom.z = safeStod(line.substr(20, 10));
        atom.element = line.substr(31, 3);
        atom.element.erase(std::remove_if(atom.element.begin(), atom.element.end(), ::isspace), atom.element.end());
        atom.inRing = false;
        atom.inAromaticRing = false;
        atom.charge = 0;
        atom.massDifference = 0;
        atom.stereoParity = 0;
        atom.hydrogenCount = 0;
        atom.stereoCount = 0;
        atom.valence = 0;
        atom.extraInfo = "";
        // 解析额外的信息
        atom.massDifference = (line.length() >= 35) ? safeStoi(line.substr(34, 2)) : 0;
        atom.charge = (line.length() >= 38) ? safeStoi(line.substr(36, 3)) : 0;
        atom.stereoParity = (line.length() >= 41) ? safeStoi(line.substr(39, 3)) : 0;
        atom.hydrogenCount = (line.length() >= 44) ? safeStoi(line.substr(42, 3)): 0; // 减1以获得实际数量
        atom.stereoCount = (line.length() >= 47) ? safeStoi(line.substr(45, 3)) : 0;
        atom.valence = (line.length() >= 50) ? safeStoi(line.substr(48, 3)) : 0;
        atom.extraInfo = (line.length() > 51) ? line.substr(51) : "";

        elementCount[atom.element]++;
        atom.name = atom.element + std::to_string(elementCount[atom.element]);
        
        mol.atoms.push_back(atom);
    }

    // Read bonds
    for (int i = 0; i < bondCount; ++i) {
        std::getline(input, line);
        SDFBond bond;
        bond.atom1 = safeStoi(line.substr(0, 3)) - 1;
        bond.atom2 = safeStoi(line.substr(3, 3)) - 1;
        bond.type = safeStoi(line.substr(6, 3));
        bond.rotatable = (bond.type == 1);
        bond.inRing = false;
        bond.isAromatic = false;
        
        // 捕获额外的信息
        if (line.length() > 9) {
            bond.exinfo = line.substr(9);
            // 移除尾部的空白字符
            bond.exinfo.erase(bond.exinfo.find_last_not_of(" \t\n\r") + 1);
        }

        mol.bonds.push_back(bond);
        mol.atoms[bond.atom1].bonds.push_back(bond.atom2);
        mol.atoms[bond.atom2].bonds.push_back(bond.atom1);
    }

    // Read M lines until M END
    std::stringstream mInfoStream;
    while (std::getline(input, line) && !input.eof()) {
        if (line == "M  END") {
            break;
        }
        mInfoStream << line << "\n";
    }
    mol.MInfo = mInfoStream.str();

    // Read property lines until $$$$ or end of file
    std::stringstream mPropStream;
    while (std::getline(input, line) && !input.eof()) {
        if (line == "$$$$") {
            break;
        }
        mPropStream << line << "\n";
    }
    mol.MProp = mPropStream.str();
    
    // 初始化环计算标志
    mol.ringsComputed = false;
    
    // 在解析后立即计算所有环
    findAllRings(mol);
    mol.ringsComputed = true;
    
    return mol;
}

// 从文件解析SDF
SDFMolecule parseSDF(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Unable to open file: " + filename);
    }
    return parseSDFImpl(file);
}

// 从字符串内容解析SDF
SDFMolecule parseSDFContent(const std::string& content) {
    std::istringstream iss(content);
    return parseSDFImpl(iss);
}

// 添加这个辅助函数来计算四面体体积
double tetrahedronVolume(const SDFAtom& a, const SDFAtom& b, const SDFAtom& c, const SDFAtom& d) {
    double ax = b.x - a.x, ay = b.y - a.y, az = b.z - a.z;
    double bx = c.x - a.x, by = c.y - a.y, bz = c.z - a.z;
    double cx = d.x - a.x, cy = d.y - a.y, cz = d.z - a.z;
    return std::abs(ax * (by * cz - bz * cy) + ay * (bz * cx - bx * cz) + az * (bx * cy - by * cx)) / 6.0;
}

bool NIsReceptor(const SDFMolecule& mol, int atomId) {
    const SDFAtom& a = mol.atoms[atomId];
    
    // 检查是否是N原子
    if (a.element != "N") return false;

    // 计算与N原子相连的原子数量
    int bondCount = a.bonds.size();

    // 规则1: 如果N原子只连接一个或两个键，则是氢键受体
    if (bondCount <= 2) return true;

    // 规则3: 如果N原子连接4个键，不是氢键受体
    if (bondCount == 4) return false;

    // 规则2: 如果N原子连接3个键
    if (bondCount == 3) {
        // 对于sp3杂化N：3个键=有孤对电子(受体)
        // 对于sp2杂化N（如在芳香环中）：需要进一步判断
        
        // 检查N是否在芳香环中
        if (a.inAromaticRing) {
            // 芳香环中的N原子（如吡啶）孤对电子在sp2轨道中，不参与π共轭，是受体
            // 但吡咯类N的孤对电子参与π共轭，不是受体
            // 检查连接的原子是否都在环中（吡咯N连接2个环原子+1个H）
            int neighborsInRing = 0;
            for (int neighborId : a.bonds) {
                if (mol.atoms[neighborId].inRing) {
                    neighborsInRing++;
                }
            }
            // 如果3个邻居都在环中，可能是季铵盐或特殊结构，不是受体
            if (neighborsInRing == 3) return false;
            // 否则，通常是受体（如吡啶）
            return true;
        }
        
        // 非芳香环中的sp3 N，3个键意味着有孤对电子
        const SDFAtom& neighbor1 = mol.atoms[a.bonds[0]];
        const SDFAtom& neighbor2 = mol.atoms[a.bonds[1]];
        const SDFAtom& neighbor3 = mol.atoms[a.bonds[2]];

        // 计算四面体体积来判断几何构型
        double volume = tetrahedronVolume(a, neighbor1, neighbor2, neighbor3);

        // 使用更合理的阈值：平面构型体积接近0
        // 四面体构型体积较大
        const double VOLUME_THRESHOLD = 0.10;

        // 如果体积大于阈值，则认为N原子是四面体构型(sp3)，有孤对电子，是氢键受体
        // 如果体积小于阈值，则是平面构型(sp2)，孤对电子在p轨道中参与共轭，可能是酰胺等
        return volume > VOLUME_THRESHOLD;
    }

    // 默认情况（不应该到达这里）
    return false;
}


struct PairHash {
    template <class T1, class T2>
    std::size_t operator () (const std::pair<T1, T2>& p) const {
        auto h1 = std::hash<T1>{}(p.first);
        auto h2 = std::hash<T2>{}(p.second);
        return h1 ^ h2;
    }
};

// BFS辅助函数：找到从start到end的最短路径（不使用排除边）
std::vector<int> bfsShortestPath(const SDFMolecule& mol, int start, int end, 
                                   int excludeAtom1, int excludeAtom2) {
    if (start == end) return {};
    
    std::vector<int> parent(mol.atoms.size(), -1);
    std::vector<bool> visited(mol.atoms.size(), false);
    std::queue<int> q;
    
    q.push(start);
    visited[start] = true;
    
    while (!q.empty()) {
        int current = q.front();
        q.pop();
        
        for (int neighbor : mol.atoms[current].bonds) {
            // 跳过排除的边
            if ((current == excludeAtom1 && neighbor == excludeAtom2) ||
                (current == excludeAtom2 && neighbor == excludeAtom1)) {
                continue;
            }
            
            if (!visited[neighbor]) {
                visited[neighbor] = true;
                parent[neighbor] = current;
                
                if (neighbor == end) {
                    // 重构路径
                    std::vector<int> path;
                    int node = end;
                    while (node != -1) {
                        path.push_back(node);
                        node = parent[node];
                    }
                    std::reverse(path.begin(), path.end());
                    return path;
                }
                
                q.push(neighbor);
            }
        }
    }
    
    return {}; // 没有找到路径
}

void findAllRings(SDFMolecule& mol) {
    // 使用基于最短路径的环检测算法
    // 对于每条边(u,v)，找到u到v的不使用该边的最短路径
    // 如果存在这样的路径，则形成一个环
    
    std::set<std::vector<int>> uniqueRings; // 使用set去重
    std::set<std::pair<int,int>> processedEdges;
    
    for (const SDFBond& bond : mol.bonds) {
        int u = bond.atom1;
        int v = bond.atom2;
        
        // 标准化边表示（确保较小索引在前）
        std::pair<int,int> edge = (u < v) ? std::make_pair(u, v) : std::make_pair(v, u);
        
        if (processedEdges.count(edge)) continue;
        processedEdges.insert(edge);
        
        // 找到从u到v的不使用(u,v)边的最短路径
        std::vector<int> path = bfsShortestPath(mol, u, v, u, v);
        
        if (!path.empty() && path.size() >= 3) {
            // 找到一个环：path + [u]
            std::vector<int> ring = path;
            
            // 标准化环表示：找到最小元素，旋转使得最小元素在前面
            auto minIt = std::min_element(ring.begin(), ring.end());
            int minPos = std::distance(ring.begin(), minIt);
            std::rotate(ring.begin(), ring.begin() + minPos, ring.end());
            
            // 检查是否需要反转（确保字典序最小）
            std::vector<int> reversedRing = ring;
            std::reverse(reversedRing.begin() + 1, reversedRing.end());
            if (reversedRing < ring) {
                ring = reversedRing;
            }
            
            uniqueRings.insert(ring);
        }
    }
    
    // 转换为向量
    mol.rings.assign(uniqueRings.begin(), uniqueRings.end());

    // Set inRing property for atoms and bonds
    for (const auto& ring : mol.rings) {
        for (int atomIndex : ring) {
            mol.atoms[atomIndex].inRing = true;
        }
        for (size_t i = 0; i < ring.size(); ++i) {
            int atom1 = ring[i];
            int atom2 = ring[(i + 1) % ring.size()];
            for (SDFBond& bond : mol.bonds) {
                if ((bond.atom1 == atom1 && bond.atom2 == atom2) ||
                    (bond.atom1 == atom2 && bond.atom2 == atom1)) {
                    bond.inRing = true;
                    break;
                }
            }
        }
    }
}

void identifyAromaticRings(SDFMolecule& mol) {
    for (auto& ring : mol.rings) {
        if (ring.size() == 5 || ring.size() == 6) {
            bool isAromatic = true;
            for (int atomId : ring) {
                const SDFAtom& atom = mol.atoms[atomId];
                if (atom.element != "C" && atom.element != "N" && atom.element != "O" && atom.element != "S") {
                    isAromatic = false;
                    break;
                }
            }
            if (isAromatic) {
                for (int atomId : ring) {
                    mol.atoms[atomId].inAromaticRing = true;
                }
                for (size_t i = 0; i < ring.size(); ++i) {
                    int atom1 = ring[i];
                    int atom2 = ring[(i + 1) % ring.size()];
                    for (SDFBond& bond : mol.bonds) {
                        if ((bond.atom1 == atom1 && bond.atom2 == atom2) ||
                            (bond.atom1 == atom2 && bond.atom2 == atom1)) {
                            bond.isAromatic = true;
                            break;
                        }
                    }
                }
            }
        }
    }
}

bool is_aromatic_bond(const SDFMolecule& mol, SDFBond& bond) {
    // 获取键连接的两个原子
    const SDFAtom& atom1 = mol.atoms[bond.atom1];
    const SDFAtom& atom2 = mol.atoms[bond.atom2];

    // 检查原子是否可能是芳香环的一部分（扩展支持更多元素）
    auto is_aromatic_element = [](const std::string& element) {
        return element == "C" || element == "N" || element == "O" || 
               element == "S" || element == "P" || element == "Se" || element == "As";
    };

    if (!is_aromatic_element(atom1.element) || !is_aromatic_element(atom2.element)) {
        return false;
    }

    // 检查是否在一个环中
    if (!bond.inRing) {
        return false;
    }

    // 找到包含这个键的环
    auto find_ring = [&mol, &bond]() -> const std::vector<int>* {
        for (const auto& ring : mol.rings) {
            for (size_t i = 0; i < ring.size(); ++i) {
                int j = (i + 1) % ring.size();
                if ((ring[i] == bond.atom1 && ring[j] == bond.atom2) ||
                    (ring[i] == bond.atom2 && ring[j] == bond.atom1)) {
                    return &ring;
                }
            }
        }
        return nullptr;
    };

    const std::vector<int>* ring = find_ring();
    if (!ring) {
        return false;
    }

    // 检查环的大小（芳香环通常是5、6或7元环）
    size_t ring_size = ring->size();
    if (ring_size < 5 || ring_size > 7) {
        return false;
    }

    // 计算环中的π电子数，改进的Hückel规则判断
    int pi_electrons = 0;
    for (size_t i = 0; i < ring->size(); ++i) {
        const SDFAtom& atom = mol.atoms[(*ring)[i]];
        
        if (atom.element == "C") {
            // sp2杂化的C贡献1个π电子
            if (atom.bonds.size() == 3) {
                pi_electrons += 1;
            } else if (atom.bonds.size() == 2) {
                // 卡宾类型，贡献2个π电子
                pi_electrons += 2;
            } else {
                return false;
            }
        } else if (atom.element == "N") {
            if (atom.bonds.size() == 2) {
                // sp2杂化的N（如吡咯型）贡献2个π电子
                pi_electrons += 2;
            } else if (atom.bonds.size() == 3) {
                // 检查是否是吡啶型N
                if (!NIsReceptor(mol, (*ring)[i])) {
                    // 吡咯型N，孤对电子参与共轭
                    pi_electrons += 2;
                } else {
                    // 吡啶型N，孤对电子不参与共轭
                    pi_electrons += 1;
                }
            } else if (atom.bonds.size() == 4) {
                // 季铵盐
                pi_electrons += 2;
            } else {
                return false;
            }
        } else if (atom.element == "O") {
            if (atom.bonds.size() == 2) {
                pi_electrons += 2;
            } else {
                return false;
            }
        } else if (atom.element == "S") {
            if (atom.bonds.size() == 2) {
                pi_electrons += 2;
            } else {
                pi_electrons += 2;
            }
        } else if (atom.element == "P") {
            pi_electrons += 2;
        } else if (atom.element == "Se" || atom.element == "As") {
            pi_electrons += 2;
        }
    }

    // 检查是否满足Hückel规则 (4n+2)
    bool aromatic = false;
    if (ring_size == 5 || ring_size == 6 || ring_size == 7) {
        aromatic = (pi_electrons == 6) || (pi_electrons == 10);
    }
    
    bond.isAromatic = (aromatic || bond.type == 4);
    return aromatic;
}

bool isAromaticCarbon(const SDFMolecule& mol, const SDFAtom& atom) {
    if (atom.element != "C") return false;
    return atom.inAromaticRing && atom.bonds.size() < 4;
}

void setAtomTypes(SDFMolecule& mol) {
    identifyAromaticRings(mol);
    for (SDFAtom& atom : mol.atoms) {     
        atom.type = atom.element.size() < 2 ? atom.element + ' ' : atom.element.substr(0, 2);
        if (atom.element == "H") {
            if (atom.bonds.size() != 1) {
                std::cerr << "H atom (id: " << atom.id + 1 << ") has " << atom.bonds.size() << " bonds." << std::endl;
                continue;
            }
            int bondIndex = atom.bonds[0];
            const SDFAtom& connectedAtom = mol.atoms[bondIndex];
            if (connectedAtom.element == "O" || connectedAtom.element == "N") {
                atom.type = "HD";
            } else {
                atom.type = "H ";
            }
        } else if (atom.element == "O") {
            bool connectedToH = false;
            for (int neighborId : atom.bonds) {
                const SDFAtom& connectedAtom = mol.atoms[neighborId];
                if (connectedAtom.element == "H") {
                    connectedToH = true;
                    break;
                }
            }
            atom.type = connectedToH ? "OA" : "OA";
        } else if (atom.element == "N") {
            bool connectedToH = false;
            for (int neighborId : atom.bonds) {
                const SDFAtom& connectedAtom = mol.atoms[neighborId];
                if (connectedAtom.element == "H") {
                    connectedToH = true;
                    break;
                }
            }
            if (NIsReceptor(mol, atom.id)) {
                atom.type = "NA";
            } else if (connectedToH) {
                atom.type = "N ";
            }
        } else if (atom.element == "S") {
            if (atom.bonds.size() <= 2) {
                atom.type = "SA";
            }
        } else if (atom.element == "P") {
            // 磷原子的类型判断
            if (atom.inRing) {
                atom.type = "PA";  // 环内磷原子
            } else {
                atom.type = "P ";  // 环外磷原子
            }
        } else if (atom.element == "F") {
            // 氟原子
            bool connectedToH = false;
            for (int neighborId : atom.bonds) {
                const SDFAtom& connectedAtom = mol.atoms[neighborId];
                if (connectedAtom.element == "H") {
                    connectedToH = true;
                    break;
                }
            }
            atom.type = connectedToH ? "FD" : "FA";
        } else if (atom.element == "Cl") {
            // 氯原子
            atom.type = "Cl";
        } else if (atom.element == "Br") {
            // 溴原子
            atom.type = "Br";
        } else if (atom.element == "I") {
            // 碘原子
            atom.type = "I ";
        } else if (atom.element == "C") {
            if (isAromaticCarbon(mol, atom)) {
                atom.type = "A ";
            } else {
                atom.type = "C "; // 或者其他适当的非芳香碳类型
            }
        } 
    }
}

bool is_amide_bond(const SDFMolecule& mol, const SDFBond& bond) {
    const SDFAtom& atom1 = mol.atoms[bond.atom1];
    const SDFAtom& atom2 = mol.atoms[bond.atom2];

    // 检查是否是C-N键
    if (!((atom1.element == "C" && atom2.element == "N") || 
          (atom1.element == "N" && atom2.element == "C"))) {
        return false;
    }

    // 确定哪个是C原子，哪个是N原子
    const SDFAtom& C_atom = (atom1.element == "C") ? atom1 : atom2;
    const SDFAtom& N_atom = (atom1.element == "N") ? atom1 : atom2;

    // 检查C原子是否与O形成双键
    bool C_double_bonded_to_O = false;
    for (int neighbor_atom_id : C_atom.bonds) {
        const SDFAtom& neighbor_atom = mol.atoms[neighbor_atom_id];
        // 找到连接这两个原子的键
        auto it = std::find_if(mol.bonds.begin(), mol.bonds.end(), 
            [&](const SDFBond& b) { 
                return (b.atom1 == C_atom.id && b.atom2 == neighbor_atom_id) || 
                       (b.atom2 == C_atom.id && b.atom1 == neighbor_atom_id); 
            });
        if (it != mol.bonds.end() && neighbor_atom.element == "O" && it->type == 2) {
            C_double_bonded_to_O = true;
            break;
        }
    }
    if (!C_double_bonded_to_O) { return false; } 

    // 检查N原子的连接情况
    bool N_bonded_to_H = false;
    for (int neighbor_atom_id : N_atom.bonds) {
        const SDFAtom& neighbor_atom = mol.atoms[neighbor_atom_id];
        // 找到连接这两个原子的键
        auto it = std::find_if(mol.bonds.begin(), mol.bonds.end(), 
            [&](const SDFBond& b) { 
                return (b.atom1 == N_atom.id && b.atom2 == neighbor_atom_id) || 
                       (b.atom2 == N_atom.id && b.atom1 == neighbor_atom_id); 
            });
        int bond_type = (it != mol.bonds.end()) ? it->type : -1;
        if (neighbor_atom.element == "H") {
            N_bonded_to_H = true;
            break;
        } 
    }

    // 判断是否为酰胺键
    // 酰胺键的条件：
    // 1. C=O双键（已在上面验证）
    // 2. C-N单键（通过bond.type == 1验证）
    // 3. N连接H原子（或至少可以形成氢键）
    // 4. N不在环中（避免将环内酰胺误判）
    bool is_amide = (C_double_bonded_to_O && bond.type == 1 && !N_atom.inRing);
    return is_amide;
}


// 检查是否是终端甲基或三氟甲基
bool isTerminalMethylOrTrifluoromethyl(const SDFMolecule& mol, int atomId) {
    const SDFAtom& atom = mol.atoms[atomId];
    if (atom.element != "C") return false;
    
    int hCount = 0;
    int fCount = 0;
    
    for (int neighborId : atom.bonds) {
        const SDFAtom& neighbor = mol.atoms[neighborId];
        if (neighbor.element == "H") hCount++;
        else if (neighbor.element == "F") fCount++;
    }
    
    return (hCount == 3) || (fCount == 3);
}

void determineRotatableBonds(SDFMolecule& mol) {
    // 环已在parseSDFImpl中计算，无需重复调用
    for (SDFBond& bond : mol.bonds) {
        if (bond.type != 1) {
            bond.rotatable = false;
            continue;
        }

        if (bond.inRing) {
            bond.rotatable = false;
            continue;
        }

        if (is_aromatic_bond(mol, bond)) {
            bond.rotatable = false;
            continue;
        }
        if (is_amide_bond(mol, bond)) {
            bond.rotatable = false;
            continue;
        }
        // Check if either atom has only one bond or is a terminal methyl/trifluoromethyl group
        if (mol.atoms[bond.atom1].bonds.size() == 1 || mol.atoms[bond.atom2].bonds.size() == 1 ||
            isTerminalMethylOrTrifluoromethyl(mol, bond.atom1) || isTerminalMethylOrTrifluoromethyl(mol, bond.atom2)) {
            bond.rotatable = false;
        }
    }
    setAtomTypes(mol);
}

std::vector<std::vector<int>> splitMolecule(const SDFMolecule& mol) {
    std::vector<std::vector<int>> fragments;
    std::vector<bool> visited(mol.atoms.size(), false);

    for (size_t i = 0; i < mol.atoms.size(); ++i) {
        if (!visited[i]) {
            std::vector<int> fragment;
            std::vector<int> stack = {static_cast<int>(i)};

            while (!stack.empty()) {
                int current = stack.back();
                stack.pop_back();

                if (!visited[current]) {
                    visited[current] = true;
                    fragment.push_back(current);

                    for (int neighbor : mol.atoms[current].bonds) {
                        if (!visited[neighbor]) {
                            bool rotatable = false;
                            for (const SDFBond& bond : mol.bonds) {
                                if ((bond.atom1 == current && bond.atom2 == neighbor) ||
                                    (bond.atom2 == current && bond.atom1 == neighbor)) {
                                    rotatable = bond.rotatable;
                                    break;
                                }
                            }
                            if (!rotatable) {
                                stack.push_back(neighbor);
                            }
                        }
                    }
                }
            }

            fragments.push_back(fragment);
        }
    }

    return fragments;
}

std::vector<int> findLargestFragmentAfterCut(const SDFMolecule& mol, const std::vector<int>& fragment) {
    std::vector<int> largestCutFragment;
    std::unordered_set<int> fragmentSet(fragment.begin(), fragment.end());
    std::unordered_set<int> allAtoms;
    for (size_t i = 0; i < mol.atoms.size(); ++i) {
        allAtoms.insert(i);
    }

    for (const auto& bond : mol.bonds) {
        if (bond.rotatable && 
            (fragmentSet.count(bond.atom1) || fragmentSet.count(bond.atom2))) {
            
            // This rotatable bond connects to the fragment
            std::vector<bool> visited(mol.atoms.size(), false);
            std::vector<int> cutFragment;

            // DFS to find connected component after cutting the bond
            std::function<void(int)> dfs = [&](int atom) {
                visited[atom] = true;
                if (!fragmentSet.count(atom)) {
                    cutFragment.push_back(atom);
                }
                for (int neighbor : mol.atoms[atom].bonds) {
                    if (!visited[neighbor] && 
                        !((atom == bond.atom1 && neighbor == bond.atom2) ||
                          (atom == bond.atom2 && neighbor == bond.atom1))) {
                        dfs(neighbor);
                    }
                }
            };

            // Start DFS from the atom not in the fragment
            int startAtom = fragmentSet.count(bond.atom1) ? bond.atom2 : bond.atom1;
            dfs(startAtom);

            if (cutFragment.size() > largestCutFragment.size()) {
                largestCutFragment = cutFragment;
            }
        }
    }

    return largestCutFragment;
}

std::vector<std::pair<int, int>> findConnectedRotatableBonds(const SDFMolecule& mol, const std::vector<int>& fragment) {
    std::vector<std::pair<int, int>> connectedBonds;
    std::unordered_set<int> fragmentSet(fragment.begin(), fragment.end());
    
    for (const auto& bond : mol.bonds) {
        if (bond.rotatable && 
            ((fragmentSet.count(bond.atom1) && !fragmentSet.count(bond.atom2)) ||
             (fragmentSet.count(bond.atom2) && !fragmentSet.count(bond.atom1)))) {
            if (bond.atom1 < 0 || bond.atom1 >= mol.atoms.size() || 
                bond.atom2 < 0 || bond.atom2 >= mol.atoms.size()) {
                std::cerr << "Error: Invalid atom indices in bond: " << bond.atom1 << " " << bond.atom2 << std::endl;
                continue;
            }
            connectedBonds.push_back({bond.atom1, bond.atom2});
        }
    }
    
    return connectedBonds;
}

std::string outputBranchStructure(const SDFMolecule& mol, const std::vector<int>& rootFragment, 
                            std::unordered_set<int>& processedAtoms) {
    std::stringstream outstream;
    std::function<void(int, int)> processBranch = [&](int parentAtom, int startAtom) {
        if (processedAtoms.count(startAtom) > 0) return;

        std::vector<int> branchFragment;
        std::function<void(int)> dfs = [&](int atom) {
            if (processedAtoms.count(atom) == 0) {
                branchFragment.push_back(atom);
                processedAtoms.insert(atom);
                for (int neighbor : mol.atoms[atom].bonds) {
                    bool isRotatable = false;
                    for (const auto& bond : mol.bonds) {
                        if ((bond.atom1 == atom && bond.atom2 == neighbor) ||
                            (bond.atom2 == atom && bond.atom1 == neighbor)) {
                            isRotatable = bond.rotatable;
                            break;
                        }
                    }
                    if (!isRotatable && processedAtoms.count(neighbor) == 0) {
                        dfs(neighbor);
                    }
                }
            }
        };

        dfs(startAtom);

        if (!branchFragment.empty()) {
            outstream << "BRANCH " << parentAtom + 1 << " " << startAtom + 1 << std::endl;

            // Output atoms in this branch
            for (int atom : branchFragment) {
                const SDFAtom& currentAtom = mol.atoms[atom];
                outstream << std::left << std::setw(6) << "ATOM"
                    << std::right << std::setw(5) << atom + 1
                    << " "
                    << std::left << std::setw(4) << currentAtom.name
                    << " "
                    << std::left << std::setw(4) << "LIG"
                    << std::left << std::setw(1) << "A"
                    << std::right << std::setw(4) << "1"
                    << "    " //A Char, iCode
                    << std::right << std::fixed << std::setprecision(3)
                    << std::setw(8) << currentAtom.x
                    << std::setw(8) << currentAtom.y
                    << std::setw(8) << currentAtom.z
                    << std::setw(6) << std::fixed << std::setprecision(2) << 1.00  // Occupancy
                    << std::setw(6) << std::fixed << std::setprecision(2) << 0.00  // Temperature factor
                    << std::setw(10) << std::fixed << std::setprecision(4) << 0.00  // Empty space
                    << " "
                    << std::left << std::setw(2) << currentAtom.type
                    << std::endl;
            }

            // Process sub-branches
            for (int atom : branchFragment) {
                for (int neighbor : mol.atoms[atom].bonds) {
                    if (processedAtoms.count(neighbor) == 0) {
                        bool isRotatable = false;
                        for (const auto& bond : mol.bonds) {
                            if ((bond.atom1 == atom && bond.atom2 == neighbor) ||
                                (bond.atom2 == atom && bond.atom1 == neighbor)) {
                                isRotatable = bond.rotatable;
                                break;
                            }
                        }
                        if (isRotatable) {
                            processBranch(atom, neighbor);
                        }
                    }
                }
            }

            outstream << "ENDBRANCH " << parentAtom + 1 << " " << startAtom + 1 << std::endl;
        }
    };

    // Process all branches connected to the root fragment
    for (int rootAtom : rootFragment) {
        for (int neighbor : mol.atoms[rootAtom].bonds) {
            if (processedAtoms.count(neighbor) == 0) {
                bool isRotatable = false;
                for (const auto& bond : mol.bonds) {
                    if ((bond.atom1 == rootAtom && bond.atom2 == neighbor) ||
                        (bond.atom2 == rootAtom && bond.atom1 == neighbor)) {
                        isRotatable = bond.rotatable;
                        break;
                    }
                }
                if (isRotatable) {
                    processBranch(rootAtom, neighbor);
                }
            }
        }
    }

     //std::cout << "Exiting outputBranchStructure" << std::endl;
     return outstream.str();
}

std::string outputOptimalFragmentAsPDB(const SDFMolecule& mol) {
    auto fragments = splitMolecule(mol);
    std::vector<std::vector<int>> optimalFragments;
    int maxConnectedRotatableBonds = -1;
    std::vector<int> largestCutFragment;

    // First rule: find fragments with the most connected rotatable bonds
    for (const auto& fragment : fragments) {
        int connectedRotatableBonds = 0;
        for (const auto& bond : mol.bonds) {
            if (bond.rotatable && 
                ((std::find(fragment.begin(), fragment.end(), bond.atom1) != fragment.end() &&
                  std::find(fragment.begin(), fragment.end(), bond.atom2) == fragment.end()) ||
                 (std::find(fragment.begin(), fragment.end(), bond.atom2) != fragment.end() &&
                  std::find(fragment.begin(), fragment.end(), bond.atom1) == fragment.end()))) {
                connectedRotatableBonds++;
            }
        }
        
        
        if (connectedRotatableBonds > maxConnectedRotatableBonds) {
            maxConnectedRotatableBonds = connectedRotatableBonds;
            optimalFragments.clear();
            optimalFragments.push_back(fragment);
        } else if (connectedRotatableBonds == maxConnectedRotatableBonds) {
            optimalFragments.push_back(fragment);
        }
    }
    std::cout << "optimalFragments.size() by rule 1= " << optimalFragments.size() << std::endl;
    // Second rule: if multiple fragments have the same max connected rotatable bonds,
    // choose the one that produces the largest cut fragment
    std::vector<int> optimalFragment;
    if (optimalFragments.size() > 1) {
        int maxCutSize = INT_MAX;
        for (const auto& fragment : optimalFragments) {
            auto cutFragment = findLargestFragmentAfterCut(mol, fragment);
            if (cutFragment.size() < maxCutSize) {
                maxCutSize = cutFragment.size();
                largestCutFragment = cutFragment;
                optimalFragment = fragment;
            }
        }
    } else if (!optimalFragments.empty()) {
        optimalFragment = optimalFragments[0];
        largestCutFragment = findLargestFragmentAfterCut(mol, optimalFragment);
    }

    if (optimalFragment.empty()) {
        std::cout << "No valid fragment found. Using the largest original fragment." << std::endl;
        optimalFragment = *std::max_element(fragments.begin(), fragments.end(),
            [](const std::vector<int>& a, const std::vector<int>& b) {
                return a.size() < b.size();
            });
    }

    std::stringstream outstream;

    outstream << "REMARK   1 SDF2PDBQT FOR " << mol.name << std::endl;
    outstream << "REMARK   2 Generated by WATVINA"   << std::endl;
    outstream << "ROOT" << std::endl;
    std::unordered_set<int> processedAtoms;

    // Output the optimal fragment (ROOT) atoms
    for (int i : optimalFragment) {
        if (i < 0 || i >= mol.atoms.size()) {
            std::cerr << "Error: Invalid atom index " << i << std::endl;
            continue;
        }
        const SDFAtom& atom = mol.atoms[i];
        outstream << std::left << std::setw(6) << "ATOM"
                << std::right << std::setw(5) << i + 1  // Use atom index + 1 as serial
                << " "
                << std::left << std::setw(4) << atom.name
                << " "
                << std::left << std::setw(4) << "LIG"
                << std::left << std::setw(1) << "A"
                << std::right << std::setw(4) << "1"
                << "    " //A Char, iCode
                << std::right << std::fixed << std::setprecision(3)
                << std::setw(8) << atom.x
                << std::setw(8) << atom.y
                << std::setw(8) << atom.z
                << std::setw(6) << std::fixed << std::setprecision(2) << 1.00  // Occupancy
                << std::setw(6) << std::fixed << std::setprecision(2) << 0.00  // Temperature factor
                << std::setw(10) << std::fixed << std::setprecision(4) << 0.00
                << " "
                << std::left << std::setw(2) << atom.type
                << std::endl;
        processedAtoms.insert(i);
    }

    outstream << "ENDROOT" << std::endl;

    outstream << outputBranchStructure(mol, optimalFragment, processedAtoms);
    int totalRotatableBonds = 0;
    for (const SDFBond& bond : mol.bonds) {
        if (bond.rotatable) {
            totalRotatableBonds++;
        }
    }    
    outstream << "TORSDOF " << totalRotatableBonds << std::endl;
    return outstream.str();
}

std::string SDF2PDBQT(const std::string& filename) {
    SDFMolecule mol = parseSDF(filename);
    // 环已在parseSDFImpl中计算，无需重复调用
    determineRotatableBonds(mol);
    setAtomTypes(mol);
    return outputOptimalFragmentAsPDB(mol);
}

std::pair<SDFMolecule, std::string> SDF2MOLPDBQT(const std::string& filename) {
    std::pair<SDFMolecule, std::string> result;
    SDFMolecule mol = parseSDF(filename);
    // 环已在parseSDFImpl中计算，无需重复调用
    determineRotatableBonds(mol);
    setAtomTypes(mol);
    std::string contents = outputOptimalFragmentAsPDB(mol);
    result.second = contents;
    result.first = mol;
    return result;
}

std::pair<SDFMolecule, std::string> SDF2MolPdbqtFromContent(const std::string& sdfContent) {
    std::pair<SDFMolecule, std::string> result;
    SDFMolecule mol = parseSDFContent(sdfContent);
    // 环已在parseSDFImpl中计算，无需重复调用
    determineRotatableBonds(mol);
    setAtomTypes(mol);
    std::string contents = outputOptimalFragmentAsPDB(mol);
    result.second = contents;
    result.first = mol;
    return result;
}


// 更新原子坐标的函数
void updateAtomCoordinates(SDFMolecule& mol, const std::map<int, std::array<double, 3>>& newCoordinates) {
    for (auto& atom : mol.atoms) {
        auto it = newCoordinates.find(atom.id);
        if (it != newCoordinates.end()) {
            atom.x = it->second[0];
            atom.y = it->second[1];
            atom.z = it->second[2];
        }
    }
}

// 将分子信息写入 SDF 格式的字符串
// 支持可选的坐标更新和备注信息
std::string writeMoleculeToSDFString(const SDFMolecule& mol, 
                                     const std::map<int, std::array<double, 3>>& newCoordinates = {},
                                     const std::string& remarkSDF = "") {
    SDFMolecule tmpmol = mol;
    updateAtomCoordinates(tmpmol, newCoordinates);
    std::stringstream ss;

    // 写入分子名称
    ss << tmpmol.name << "\n";
    ss << "  Generated by WATVINA\n";
    ss << "\n";

    // 写入原子数和键数
    ss << std::setw(3) << tmpmol.atoms.size() << std::setw(3) << tmpmol.bonds.size() << "  0  0  0  0  0  0  0  0999 V2000\n";

    // 写入原子信息
    for (const auto& atom : tmpmol.atoms) {
        ss << std::fixed << std::setprecision(4)
           << std::setw(10) << atom.x
           << std::setw(10) << atom.y
           << std::setw(10) << atom.z
           << " "<< std::left << std::setw(3) << atom.element
           << std::right << std::setw(2) << atom.massDifference
           << std::right << std::setw(3) << atom.charge
           << std::setw(3) << atom.stereoParity
           << std::setw(3) << atom.hydrogenCount
           << std::setw(3) << atom.stereoCount
           << std::setw(3) << atom.valence
           << atom.extraInfo 
           << "\n";
    }

    // 写入键信息
    for (const auto& bond : tmpmol.bonds) {
        ss << std::setw(3) << bond.atom1 + 1
           << std::setw(3) << bond.atom2 + 1
           << std::setw(3) << bond.type
           << bond.exinfo
           << "\n";
    }

    // 写入 M 信息
    ss << tmpmol.MInfo;

    // 写入结束标记
    ss << "M  END\n";

    // 写入属性信息
    ss << remarkSDF;
    ss << tmpmol.MProp;

    // 写入文件结束标记
    ss << "$$$$\n";

    return ss.str();
}


int main(int argc, char* argv[]) {
    std::string inputFile;
    std::string outputFile = "output.pdbqt";  // Default output file name

    int opt;
    while ((opt = getopt(argc, argv, "i:o:")) != -1) {
        switch (opt) {
            case 'i':
                inputFile = optarg;
                break;
            case 'o':
                outputFile = optarg;
                break;
            default:
                std::cerr << "Usage: " << argv[0] << " -i <input_file> [-o <output_file>]" << std::endl;
                return 1;
        }
    }

    if (inputFile.empty()) {
        std::cerr << "Input file is required. Use -i <input_file>" << std::endl;
        return 1;
    }

    try {
        std::cout << "Parsing SDF file: " << inputFile << std::endl;
        SDFMolecule mol = parseSDF(inputFile);
        std::cout << "Molecule parsed. Number of atoms: " << mol.atoms.size() << ", Number of bonds: " << mol.bonds.size() << std::endl;
        // 环已在parseSDFImpl中计算，无需重复调用
        std::cout << "Determining rotatable bonds..." << std::endl;
        determineRotatableBonds(mol);

        int rotatableBondCount = std::count_if(mol.bonds.begin(), mol.bonds.end(), [](const SDFBond& b) { return b.rotatable; });
        setAtomTypes(mol);
        
        std::cout << "Number of rotatable bonds: " << rotatableBondCount << std::endl;

        std::cout << "Splitting molecule into fragments..." << std::endl;
        std::vector<std::vector<int>> frags = splitMolecule(mol);
       
        std::cout << "Number of fragments: " << frags.size() << std::endl;
       
        std::string pdbqtcontent = outputOptimalFragmentAsPDB(mol);
        std::ofstream outfile(outputFile);
        if (!outfile) {
            std::cerr << "Unable to open file: " << outputFile << std::endl;
            return 1;
        }        
        outfile << pdbqtcontent;
        outfile.close();
        std::cout << "PDB file written successfully." << std::endl;

        return 0;
    } catch (const std::exception& e) {
        std::cerr << "An error occurred: " << e.what() << std::endl;
        return 1;
    } catch (...) {
        std::cerr << "An unknown error occurred." << std::endl;
        return 1;
    }
}

