#include <iostream>
#include <sstream>
#include <limits>
#include <vector>
#include <string>
#include <cmath>
#include <stack>

class Polynom
{
private:
    
    int m = 0, n = 0;
    
    std::vector<double> a;
    std::vector<double> b;
    
public:
    
    Polynom(const int new_m, const int new_n)
    {
        m = new_m;
        n = new_n;
    }
    Polynom()
    {
    }
    
    void pushBackA(const double next_coeff)
    {
        a.push_back(next_coeff);
    }
    
    void pushBackB(const double next_coeff)
    {
        b.push_back(next_coeff);
    }
    
    void pushA(const int pos, const double coeff)
    {
        a[pos] = coeff;
    }
    
    void pushB(const int pos, const double coeff)
    {
        b[pos] = coeff;
    }
    
    void copyToA(const std::vector<double> &for_copy)
    {
        a = for_copy;
    }
    
    void copyToB(const std::vector<double> &for_copy)
    {
        b = for_copy;
    }
    
    double getA(const int pos) const
    {
        return a[pos];
    }
    
    double getB(const int pos) const
    {
        return b[pos];
    }
    
    int getM() const
    {
        return m;
    }
    
    int getN() const
    {
        return n;
    }
    
    ~Polynom()
    {
        m = 0;
        n = 0;
        
        a.clear();
        b.clear();
    }
};

bool is_expression(std::stringstream &ss)
{
  // Read std::cin to stringstream
  while (true)
  {
    std::string line("");
    getline(std::cin, line);
    if(std::cin.eof())
      break;
    ss << line << " ";
  }

  // If we find '/' or '*' symbol, that mean
  // there is expression.
  if( ss.str().find('/') != std::string::npos ||
      ss.str().find('*') != std::string::npos )
  {
    return true;
  }
  return false;
}

// Is our string is one of operators.
bool is_operator(const std::string& str)
{
  return str == "/" || str == "*" || str == "+" || str == "-";
}

// Is string starts with x^
bool is_x_n_term(const std::string& str)
{
  return str.find("x^") == 0;
}

// Check that str is double
bool is_digits(const char* str)
{
    char* endptr = 0;
    strtod(str, &endptr);

    if(*endptr != '\0' || endptr == str)
        return false;
    return true;
}

// Parse expression and extract polynom coefficients
// (example: "((1/(x^1+1))*(2*x^2))" --> A={2, 0, 0 }, B = {1, 1} )
Polynom parse_expression(std::stringstream & ss)
{
  // size of vector with coefficient 
  int n = 0, m = 0;
  std::vector<double> a;
  std::vector<double> b;
  
  char ch;
  std::stringstream result;
  std::stack<std::string> st; 
  int lenght = 0;
  int lenght_after_div = -1;
  bool is_denominator = false;
  bool meet_denominator = false;
  while(ss >> ch)
  {
    switch(ch)
    {
      case '(':
      {
        ++lenght;
        if(!result.str().empty()) 
        {
          st.push(result.str());
          result.str(std::string(""));
          st.push(std::string("("));
        }
        break;
      }
      case ')':
      { 
        if(!result.str().empty())
        { 
          if(is_x_n_term(result.str()))
          {
            int pow = std::stoi(result.str().substr(2));
            std::vector<double> &link = is_denominator ? b : a;
            if(link.size() == 0)
              link.resize(pow + 1, 0);
            if(st.size() != 0 && st.top() == "*")
            {
              st.pop();
              if(st.size() != 0 && is_digits(st.top().c_str()))
              {
                double number = std::stod(st.top());
                st.pop();
                if(st.size() != 0 && st.top() == "-")
                  number *= -1;
                link[pow] = number;
              } else {
                link[pow] = 1;
              }
            } else {
              link[pow] = 1;
            }
          } 
          else if(is_digits(result.str().c_str()))
          {
            int number = std::stod(result.str());
            if(st.size() != 0 && st.top() == "-")
              number *= -1;
            std::vector<double> &link = is_denominator ? b : a;
            if(link.size() == 0)
              link.resize(1);
            link[0] = number;
          }
          st.push(result.str());
	  result.str(std::string(""));
        }
        // Meet x^ at begin.
        if( st.size() != 0 && is_x_n_term(st.top()))
        {
          int pow = std::stoi(st.top().substr(2));
          if(is_denominator)
          {
            if(b.size() == 0)
              b.resize(pow, 0);
            b[pow - 1] = pow;
          }
        }
        if(lenght == lenght_after_div)
        {
          st.push("=");
          lenght_after_div = -1;
          is_denominator = false;
          meet_denominator = true;
        }
        --lenght;
        st.push(std::string(")"));
        break;
      }
      case '/':
      case '+':
      case '-':
      case '*':
      {
        if(ch == '*' && meet_denominator)
        {
          a.clear();
        }
        if(ch == '+' || ch == '-' || ch == '/')
        {
          if(is_x_n_term(result.str()))
          {
            int pow = std::stoi(result.str().substr(2));
            std::vector<double> &link = is_denominator ? b : a;
            if(link.size() == 0)
              link.resize(pow + 1, 0);
            if(st.size() != 0 && st.top() == "*")
            {
              st.pop();
              double number = std::stod(st.top());
              st.pop();
              if(st.top() == "-")
                number *= -1;
              link[pow] = number;
            } else {
              link[pow] = 1.0;
            }
          } else if(is_digits(result.str().c_str()))
          {
            int number = std::stod(result.str());
            if(st.size() != 0 && st.top() == "-")
              number *= -1;
            std::vector<double> &link = is_denominator ? b : a;
            if(link.size() == 0)
              link.resize(1);
            link[0] = number;
          }
        }
        
        if(!result.str().empty())
        {
          st.push(result.str());
          // Meet x^ at begin.
          if(st.top().find("x^") == 0)
          { 
            std::string pow = st.top().substr(2);
          }
        }
        result.str(std::string(""));
        if(ch == '/')
        {
          lenght_after_div = lenght;
          is_denominator = true;
        }
        std::string tmp("");
        tmp += ch;
        st.push(tmp);
        break;
      }
      default:
      { 
        result << ch;
        break;
      } 
    }
    if(ss.eof())
      break;
  }
   
  Polynom p(a.size()-1, b.size()-1);
  for(size_t i = 0; i < a.size(); ++i)
  {
    p.pushBackA(a[i]);
  }
  for(size_t i = 0; i < b.size(); ++i)
  {
    p.pushBackB(b[i]);
  }
  return p;
}

double valueat(const Polynom &p, const double point)
{
    double val_a_in_point = 0;
    double val_b_in_point = 0;
    
    for(int i = 0; i <= p.getM(); ++i)
        val_a_in_point += p.getA(i) * pow(point, i);
    
    for(int i = 0; i <= p.getN(); ++i)
        val_b_in_point += p.getB(i) * pow(point, i);
    
    return (val_a_in_point / val_b_in_point);
}

Polynom differentiate(const Polynom &p)
{
    int m = 0;
    if(((p.getM() - 1) + p.getN()) > ((p.getN() - 1) + p.getM()))
        m = (p.getM() - 1) + p.getN();
    else
        m = (p.getN() - 1) + p.getM();
    
    Polynom derivative(m, p.getN() * 2);
    
    std::vector<double> der_a(m + 1, 0);
    std::vector<double> der_b(p.getN() * 2 + 1, 0);
    
    std::vector<double> der_u(p.getM());
    std::vector<double> der_v(p.getN());
    
    for(int i = 1; i <= p.getM(); ++i)
        der_u[i - 1] = p.getA(i) * i;
    
    for(int i = 1; i <= p.getN(); ++i)
        der_v[i - 1] = p.getB(i) * i;
    
    // numerator
    for(int i = 0; i < der_u.size(); ++i)
        for(int j = 0; j <= p.getN(); ++j)
            der_a[i + j] += der_u[i] * p.getB(j);        //derivative.pushA(i + j, der_u[i] * p.getB(j));         
        
    for(int i = 0; i < der_v.size(); ++i)
        for(int j = 0; j <= p.getM(); ++j)
            der_a[i + j] -= der_v[i] * p.getA(j);   //derivative.pushA(i + j, -1 * der_v[i] * p.getA(j));
    
    derivative.copyToA(der_a);  
        
    // denominator
    for(int i = 0; i <= p.getN(); ++i)
        for(int j = 0; j <= p.getN(); ++j)
            der_b[i + j] += p.getB(i) * p.getB(j);     //derivative.pushB(i + j, p.getB(i) * p.getB(j));
    
    derivative.copyToB(der_b);    
    
    return derivative;    
}

double newton(const Polynom &p, double x_start)
{
    double eps = 0.0000000001;
    double x_last = 0;
    double x_next = x_start;
    
    Polynom derivative_p = differentiate(p);
    
    do {
        x_last = x_next;
        x_next = x_last - (valueat(p, x_last) / valueat(derivative_p, x_last));
    } while(fabs(x_next - x_last) > std::numeric_limits<double>::epsilon());
    
    return x_last;
}

double root(const Polynom &p)
{
    double x0 = 1.0;
    return newton(p, x0);
}

Polynom read_data(const std::stringstream& ss)
{
  std::cin.rdbuf(ss.rdbuf());
  int m = 0, n = 0;
  std::cin >> m >> n;

  Polynom p(m, n);
  
  // Read a0, a1, ...
  for(size_t i = 0; i <= m; ++i)
  {
    double coeff = 0;
    std::cin >> coeff;
    p.pushBackA(coeff);
  }

  // Read b0, b1, ...
  for(size_t i = 0; i <= n; ++i)
  {
    double coeff = 0;
    std::cin >> coeff;
    p.pushBackB(coeff);
  }
  return p;
}

int main()
{
    std::stringstream ss("");
    
    // Input it's expression or array of coefficients?
    bool is_expr = is_expression(ss);
    
    Polynom p;
    
    // if input is expression, parse it.
    if(is_expr)
    {
      p = parse_expression(ss);
    } else 
    {
      p = read_data(ss);
    }

    // Compute x1 with Newton method for f(x)
    double x1 = root(p);
    Polynom p2 = differentiate(p);
    
    // Compute x2 with Newton method for f'(x)
    double x2 = root(p2);
    
    // Detect maximum and minimum
    Polynom p3 = differentiate(p2);
    double x3 = valueat(p3, x1);
    double x4 = valueat(p3, x2);
    
    
    if(x3 < 0)
    {
      std::cout << "x1 " << x3 << " is maximum" << std::endl;
    } else if (x3 > 0) 
    {
      std::cout << "x1 " << x3 << " is minimum" << std::endl;
    }
    
    if(x4 < 0 and x3 != x4)
    {
      std::cout << "x2 " << x4 << " is maximum" << std::endl;
    } else if (x4 > 0 && x3 != x4)
    {
      std::cout << "x2 " << x4 << " is minimum" << std::endl;
    }
    return 0;
}
