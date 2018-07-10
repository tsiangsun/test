#include <iostream>
using namespace std;

inline int positive_modulo(int i, int n);

int main (int argc, char *argv[]) {
    int a = 5;
    int b = 10;
    int c = -3;
    
    cout << "5 % 10 = " << a%b << endl;
    cout << "-3 % 10 = " << c%b << endl;
    cout << "positive_modulo(-3, 10) = " << positive_modulo(c,b) << endl;
    
    return (0);
}

inline int positive_modulo(int i, int n) {
    return (i % n + n) % n;
}
