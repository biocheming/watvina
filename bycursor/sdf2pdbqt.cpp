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

#include <stdexcept>
#include <getopt.h>

struct Atom {
    int id;
    double x, y, z;
    std::string element;
    std::vector<int> bonds;
    bool inRing;
     Atom() : x(0), y(0), z(0), inRing(false) {}  // 默认构造函数
};

struct Bond {
    int atom1, atom2;
    int type;
    bool rotatable;
    bool inRing;
    Bond() : atom1(-1), atom2(-1), type(0), rotatable(false), inRing(false) {}  // 默认构造函数

};

struct Molecule {
    std::vector<Atom> atoms;
    std::vector<Bond> bonds;
    std::vector<std::vector<int>> rings;
};

Molecule parseSDF(const std::string& filename) {
    Molecule mol;
    std::ifstream file(filename);
    std::string line;

    // Skip header
    for (int i = 0; i < 3; ++i) {
        std::getline(file, line);
    }

    // Read atom and bond counts
    std::getline(file, line);
    int atomCount = std::stoi(line.substr(0, 3));
    int bondCount = std::stoi(line.substr(3, 3));

    // Read atoms
    for (int i = 0; i < atomCount; ++i) {
        std::getline(file, line);
        Atom atom;
        atom.id = i;
        atom.x = std::stod(line.substr(0, 10));
        atom.y = std::stod(line.substr(10, 10));
        atom.z = std::stod(line.substr(20, 10));
        atom.element = line.substr(31, 2);
        atom.element.erase(atom.element.find_last_not_of(" ") + 1);
        atom.inRing = false;
        mol.atoms.push_back(atom);
    }

    // Read bonds
    for (int i = 0; i < bondCount; ++i) {
        std::getline(file, line);
        Bond bond;
        bond.atom1 = std::stoi(line.substr(0, 3)) - 1;
        bond.atom2 = std::stoi(line.substr(3, 3)) - 1;
        bond.type = std::stoi(line.substr(6, 3));
        bond.rotatable = (bond.type == 1);  // Initially set single bonds as rotatable
        bond.inRing = false;
        mol.bonds.push_back(bond);
        mol.atoms[bond.atom1].bonds.push_back(bond.atom2);
        mol.atoms[bond.atom2].bonds.push_back(bond.atom1);
    }

    return mol;
}

// 计算点到平面的距离
double distanceToPlane(const std::array<double, 3>& point, 
                       const std::array<double, 3>& planePoint1,
                       const std::array<double, 3>& planePoint2,
                       const std::array<double, 3>& planePoint3) {
    // 计算平面法向量
    double nx = (planePoint2[1] - planePoint1[1]) * (planePoint3[2] - planePoint1[2]) - 
                (planePoint2[2] - planePoint1[2]) * (planePoint3[1] - planePoint1[1]);
    double ny = (planePoint2[2] - planePoint1[2]) * (planePoint3[0] - planePoint1[0]) - 
                (planePoint2[0] - planePoint1[0]) * (planePoint3[2] - planePoint1[2]);
    double nz = (planePoint2[0] - planePoint1[0]) * (planePoint3[1] - planePoint1[1]) - 
                (planePoint2[1] - planePoint1[1]) * (planePoint3[0] - planePoint1[0]);
    
    // 归一化法向量
    double length = std::sqrt(nx*nx + ny*ny + nz*nz);
    nx /= length; ny /= length; nz /= length;

    // 计算点到平面的距离
    return std::abs(nx * (point[0] - planePoint1[0]) + 
                    ny * (point[1] - planePoint1[1]) + 
                    nz * (point[2] - planePoint1[2]));
}

bool NIsReceptor(const Molecule& mol, int atomId) {
    const Atom& a = mol.atoms[atomId];
    
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
        std::array<double, 3> nPoint = {a.x, a.y, a.z};
        std::array<std::array<double, 3>, 3> neighborPoints;

        for (int i = 0; i < 3; ++i) {
            const Atom& neighbor = mol.atoms[a.bonds[i]];
            neighborPoints[i] = {neighbor.x, neighbor.y, neighbor.z};
        }

        double distance = distanceToPlane(nPoint, neighborPoints[0], neighborPoints[1], neighborPoints[2]);

        // 如果距离大于0.3，则是氢键受体
        return distance > 0.3;
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

void findRingsDFS(const Molecule& mol, int start, int current, std::vector<bool>& visited, 
                  std::vector<int>& path, std::vector<std::vector<int>>& rings) {
    visited[current] = true;
    path.push_back(current);

    for (int neighbor : mol.atoms[current].bonds) {
        if (neighbor == start && path.size() > 2) {
            rings.push_back(path);
        } else if (!visited[neighbor]) {
            findRingsDFS(mol, start, neighbor, visited, path, rings);
        }
    }

    visited[current] = false;
    path.pop_back();
}

void findAllRings(Molecule& mol) {
    std::vector<bool> visited(mol.atoms.size(), false);
    std::vector<int> path;
    for (size_t i = 0; i < mol.atoms.size(); ++i) {
        findRingsDFS(mol, i, i, visited, path, mol.rings);
    }

    // Set inRing property for atoms and bonds
    for (const auto& ring : mol.rings) {
        for (int atomIndex : ring) {
            mol.atoms[atomIndex].inRing = true;
        }
        for (size_t i = 0; i < ring.size(); ++i) {
            int atom1 = ring[i];
            int atom2 = ring[(i + 1) % ring.size()];
            for (Bond& bond : mol.bonds) {
                if ((bond.atom1 == atom1 && bond.atom2 == atom2) ||
                    (bond.atom1 == atom2 && bond.atom2 == atom1)) {
                    bond.inRing = true;
                    break;
                }
            }
        }
    }
}

void determineRotatableBonds(Molecule& mol) {
    findAllRings(mol);

    for (Bond& bond : mol.bonds) {
        if (bond.type != 1) {
            bond.rotatable = false;
            continue;
        }

        if (bond.inRing) {
            bond.rotatable = false;
            continue;
        }

        // Check if either atom has only one bond
        if (mol.atoms[bond.atom1].bonds.size() == 1 || mol.atoms[bond.atom2].bonds.size() == 1) {
            bond.rotatable = false;
        }
    }
}

std::vector<std::vector<int>> splitMolecule(const Molecule& mol) {
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
                            for (const Bond& bond : mol.bonds) {
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

void outputLargestFragmentAsPDB(const Molecule& mol, const std::string& filename) {
    auto fragments = splitMolecule(mol);
    auto largestFragment = *std::max_element(fragments.begin(), fragments.end(),
        [](const std::vector<int>& a, const std::vector<int>& b) {
            return a.size() < b.size();
        });

    std::ofstream outfile(filename);
    for (int i : largestFragment) {
        const Atom& atom = mol.atoms[i];
        outfile << "ATOM  " 
		<< std::setw(5) << atom.id << " "
                << std::left << std::setw(4) << atom.element
                << " LIG A   1   "
                << std::right << std::fixed << std::setprecision(3) << std::setw(8) << atom.x
                << std::right << std::fixed << std::setprecision(3) << std::setw(8) << atom.y
                << std::right << std::fixed << std::setprecision(3) << std::setw(8) << atom.z
                << "  1.00  0.00           " << atom.element << std::endl;
    }
    outfile << "END" << std::endl;
}

std::vector<int> findLargestFragmentAfterCut(const Molecule& mol, const std::vector<int>& fragment) {
    std::vector<int> largestCutFragment;
    std::unordered_set<int> fragmentSet(fragment.begin(), fragment.end());
    std::unordered_set<int> allAtoms;
    for (size_t i = 0; i < mol.atoms.size(); ++i) {
        allAtoms.insert(i);
    }

    std::cout << "Analyzing fragment of size: " << fragment.size() << std::endl;

    for (const auto& bond : mol.bonds) {
        if (bond.rotatable && 
            (fragmentSet.count(bond.atom1) || fragmentSet.count(bond.atom2))) {
            std::cout << "Found rotatable bond: " << bond.atom1 << " - " << bond.atom2 << std::endl;
            
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

            std::cout << "Cut fragment size: " << cutFragment.size() << std::endl;

            if (cutFragment.size() > largestCutFragment.size()) {
                largestCutFragment = cutFragment;
            }
        }
    }

    std::cout << "Largest cut fragment size: " << largestCutFragment.size() << std::endl;

    return largestCutFragment;
}

std::vector<std::pair<int, int>> findConnectedRotatableBonds(const Molecule& mol, const std::vector<int>& fragment) {
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

void outputBranchStructure(const Molecule& mol, const std::vector<int>& rootFragment, 
                           std::ofstream& outfile, std::unordered_set<int>& processedAtoms) {
    std::cout << "Entering outputBranchStructure with root fragment size: " << rootFragment.size() << std::endl;

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
            outfile << "BRANCH " << parentAtom + 1 << " " << startAtom + 1 << std::endl;

            // Output atoms in this branch
            for (int atom : branchFragment) {
                const Atom& currentAtom = mol.atoms[atom];
                outfile << std::left << std::setw(6) << "ATOM"
                        << std::right << std::setw(5) << atom + 1
                        << "  "
                        << std::left << std::setw(4) << currentAtom.element
                        << std::left << std::setw(3) << "LIG"
                        << " "
                        << std::left << std::setw(1) << "A"
                        << std::right << std::setw(4) << "1"
                        << "    "
                        << std::right << std::fixed << std::setprecision(3)
                        << std::setw(8) << currentAtom.x
                        << std::setw(8) << currentAtom.y
                        << std::setw(8) << currentAtom.z
                        << std::setw(6) << std::fixed << std::setprecision(2) << 1.00
                        << std::setw(6) << std::fixed << std::setprecision(2) << 0.00
                        << std::string(10, ' ')
                        << std::right << std::setw(2) << currentAtom.element
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

            outfile << "ENDBRANCH " << parentAtom + 1 << " " << startAtom + 1 << std::endl;
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

    std::cout << "Exiting outputBranchStructure" << std::endl;
}

void outputOptimalFragmentAsPDB(const Molecule& mol, const std::string& filename) {
    auto fragments = splitMolecule(mol);
    std::vector<std::vector<int>> optimalFragments;
    int maxConnectedRotatableBonds = -1;
    std::vector<int> largestCutFragment;

    std::cout << "Number of fragments: " << fragments.size() << std::endl;

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
        
        std::cout << "Fragment size: " << fragment.size() 
                  << ", Connected rotatable bonds: " << connectedRotatableBonds << std::endl;
        
        if (connectedRotatableBonds > maxConnectedRotatableBonds) {
            maxConnectedRotatableBonds = connectedRotatableBonds;
            optimalFragments.clear();
            optimalFragments.push_back(fragment);
        } else if (connectedRotatableBonds == maxConnectedRotatableBonds) {
            optimalFragments.push_back(fragment);
        }
    }

    std::cout << "Fragments with max connected rotatable bonds: " << optimalFragments.size() << std::endl;

    // Second rule: if multiple fragments have the same max connected rotatable bonds,
    // choose the one that produces the largest cut fragment
    std::vector<int> optimalFragment;
    if (optimalFragments.size() > 1) {
        int maxCutSize = -1;
        for (const auto& fragment : optimalFragments) {
            auto cutFragment = findLargestFragmentAfterCut(mol, fragment);
            std::cout << "Fragment size: " << fragment.size() 
                      << ", Largest cut fragment size: " << cutFragment.size() << std::endl;
            if (cutFragment.size() > maxCutSize) {
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

    std::cout << "Optimal fragment size: " << optimalFragment.size() 
              << ", Connected rotatable bonds: " << maxConnectedRotatableBonds
              << ", Largest cut fragment size: " << largestCutFragment.size() << std::endl;

     std::ofstream outfile(filename);
    if (!outfile) {
        std::cerr << "Unable to open file: " << filename << std::endl;
        return;
    }

    std::cout << "Writing ROOT" << std::endl;
    outfile << "ROOT" << std::endl;
    std::unordered_set<int> processedAtoms;

    // Output the optimal fragment (ROOT) atoms
    for (int i : optimalFragment) {
        if (i < 0 || i >= mol.atoms.size()) {
            std::cerr << "Error: Invalid atom index " << i << std::endl;
            continue;
        }
        const Atom& atom = mol.atoms[i];
        outfile << std::left << std::setw(6) << "ATOM"
                << std::right << std::setw(5) << i + 1  // Use atom index + 1 as serial
                << "  "
                << std::left << std::setw(4) << atom.element
                << std::left << std::setw(3) << "LIG"
                << " "
                << std::left << std::setw(1) << "A"
                << std::right << std::setw(4) << "1"
                << "    "
                << std::right << std::fixed << std::setprecision(3)
                << std::setw(8) << atom.x
                << std::setw(8) << atom.y
                << std::setw(8) << atom.z
                << std::setw(6) << std::fixed << std::setprecision(2) << 1.00  // Occupancy
                << std::setw(6) << std::fixed << std::setprecision(2) << 0.00  // Temperature factor
                << std::string(10, ' ')  // Empty space
                << std::right << std::setw(2) << atom.element
                << std::endl;
        processedAtoms.insert(i);
    }

    outfile << "ENDROOT" << std::endl;

    std::cout << "Calling outputBranchStructure" << std::endl;
    outputBranchStructure(mol, optimalFragment, outfile, processedAtoms);
    std::cout << "outputBranchStructure completed" << std::endl;
    
    outfile << "TORSDOF " << maxConnectedRotatableBonds << std::endl;
    outfile << "END" << std::endl;
    std::cout << "PDB file writing completed" << std::endl;
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
        Molecule mol = parseSDF(inputFile);
        std::cout << "Molecule parsed. Number of atoms: " << mol.atoms.size() << ", Number of bonds: " << mol.bonds.size() << std::endl;

        std::cout << "Determining rotatable bonds..." << std::endl;
        determineRotatableBonds(mol);
        int rotatableBondCount = std::count_if(mol.bonds.begin(), mol.bonds.end(), [](const Bond& b) { return b.rotatable; });
        std::cout << "Number of rotatable bonds: " << rotatableBondCount << std::endl;

        std::cout << "Splitting molecule into fragments..." << std::endl;
        std::vector<std::vector<int>> frags = splitMolecule(mol);
        std::cout << "Number of fragments: " << frags.size() << std::endl;
        for (int i = 0; i < frags.size(); ++i) {
            std::cout << "Fragment " << i + 1 << " size: " << frags[i].size() << std::endl;
            for (int j = 0; j < frags[i].size(); ++j) {
                std::cout << "Atom " << j + 1 << " in fragment " << i + 1 << ": " << frags[i][j] << std::endl;
            }
        }

        std::cout << "Outputting optimal fragment as PDB to: " << outputFile << std::endl;
        outputOptimalFragmentAsPDB(mol, outputFile);
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
