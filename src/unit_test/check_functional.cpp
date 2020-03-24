#include <functional>
#include <iostream>

using namespace std;

class test
{
private:
    double data;
public:
    test();
    ~test();

    std::function<double()> func;

    void SetData(double tmp);
};

test::test()
{
    data = 1.0;
    func = [&](){return data;};
}

test::~test()
{
}

void test::SetData(double tmp)
{
    data = tmp;
}

double Callfunc(std::function<double()> _func)
{
    return _func();
}

int main(int argc, char const *argv[])
{
    test mod;
    mod.SetData(3.3);
    double data = Callfunc(mod.func);
    cout<<"data: "<<data<<endl;
    return 0;
}
