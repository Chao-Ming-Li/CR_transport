#include <iostream>
class test_class{
public:
    test_class(double &x);
private:
    double x_;
};

test_class::test_class (double &x) : x_(1.0) {
    if (x) x = x + 1.0;
    else   x = x_ + 1.0;
}

int main() {
    double x = 0.0;
    test_class obj(x);
    std::cout << "x after object creation: " << x << std::endl; // This will print 5.0
    return 0;
}