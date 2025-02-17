#include <matrix.h>

class Matrix{
    std::vector<std::vector> data;
    int rows, cols;
    public:
    Matrix(int r, int c) : rows(r), cols(c), data(r, std::vector<double>(c, 0.0)) {}

    void fill(){
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                std::cout << " i = " << i << " / j = " << std::endl;
                std::cin >> data[i][j]; 
                }
        }
    }

    void print(){
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                std::cout << data[i][j] << std::endl;
            }
        }
    }

    void set(double element, int r, int c) {
        if (r > rows || c > cols) { std::cout << "vishel za diapazon!" << std::endl; 
        break;
        }
        data[r][c] = element;
    }

    void printElement(int r, int c) {
        if (r > rows || c > cols) { std::cout << "vishel za diapazon!" << std::endl;
        break;
        }
        std::cout << data[r][c] << std::endl;
    }

    std::vector<double> gaussJordamSolve(std::vector<double> b) {
        if (b.size() == rows){
            std::vector<std::vector<double>> extended;
            for (int i = 0; i < rows; i++) {
                for (int j = 0; j < cols; j++) {
                    extended(i, std::vector<double>(j, data[i][j]))
                }
                extended[i].push_back(b[i]);
            }
            
        }
    }
};
