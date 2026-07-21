#include <iostream>
#include <fstream>
#include <iomanip>
#include <limits>
#include <algorithm>
#include <vector>
#include <array>
#include <string>
#include <filesystem>

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " filename [padding]" << std::endl;
        return 1;
    }

    // Check if the provided path is a regular file
    if (!std::filesystem::is_regular_file(argv[1])) {
        throw std::runtime_error("The provided path is not a regular file.");
    }

    std::ifstream lig(argv[1]);
    if (!lig.is_open()) {
        throw std::runtime_error("Cannot open file: " + std::string(argv[1]));
    }

    float exten = 5.0f;
    if (argc > 2) {
        exten = std::atof(argv[2]);
    }

    const size_t d = 3;
    std::array<size_t, d> positions = {30, 38, 46};
    std::array<char, d> coordinates = {'x', 'y', 'z'};
    std::vector<float> min_values(d, std::numeric_limits<float>::max());
    std::vector<float> max_values(d, std::numeric_limits<float>::lowest());

    std::string line;
    while (std::getline(lig, line)) {
        if (line.size() < 4 || (line.substr(0, 4) != "ATOM" && line.substr(0, 6) != "HETATM")) {
            continue;
        }
        for (size_t i = 0; i < d; ++i) {
            float value = std::stof(line.substr(positions[i], 8));
            min_values[i] = std::min(min_values[i], value);
            max_values[i] = std::max(max_values[i], value);
        }
    }

    std::cout << std::fixed << std::setprecision(3);
    for (size_t i = 0; i < d; ++i) {
        std::cout << "LEDOCK_" << coordinates[i] << ": "
                  << min_values[i] - exten << " "
                  << max_values[i] + exten << std::endl;
    }

    for (size_t i = 0; i < d; ++i) {
        float center = (max_values[i] + min_values[i]) * 0.5f;
        std::cout << "VINA: center_" << coordinates[i] << '=' << center << std::endl;
    }

    for (size_t i = 0; i < d; ++i) {
        float size = (max_values[i] - min_values[i]) + 2 * exten;
        std::cout << "VINA: size_" << coordinates[i] << '=' << size << std::endl;
    }

    return 0;
}