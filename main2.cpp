#include <iostream>
#include <vector>
#include <iomanip>

using namespace std;

int main() {
    int n;
    cin >> n;

    vector<vector<long long>> matrix(n, vector<long long>(n));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            cin >> matrix[i][j];
        }
    }


     for (int i = 0; i < n; ++i) {
       int row = i;
       int col = n - 1 - i;
       int sum_indices = row + col;
       long long sum = 0;
       int count = 0;
       int replace_row = -1;
       int replace_col = -1;
        if(sum_indices > 0 && sum_indices < 2*n - 2){
            for (int r = 0; r < n; ++r) {
                for (int c = 0; c < n; ++c) {
                    if (r < c && r + c > n - 1) {
                        sum += matrix[r][c];
                        count++;
                        replace_row = r;
                        replace_col = c;
                    }
                }
            }
            if (count > 0) {
                 matrix[i][n - 1 - i] = (replace_row + replace_col);
              }else{
                 matrix[i][n - 1 - i] = -1;
            }
        }else{
             matrix[i][n-1-i] = -1;
        }
    }


    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            cout << matrix[i][j] << (j == n - 1 ? "" : " ");
        }
        cout << endl;
    }

    return 0;
}