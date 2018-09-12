#include<bits/stdc++.h>
using namespace std;
const int MAX = 1e7;
const int C = 1e6;
 
#define x first.first
#define fx first.second
#define dif second
 
bitset<MAX> sign;
vector<pair<double, double>> values;
int iteracoes;
double EPS = 1e-5;
int contador, casas;
 
vector< pair< pair<double,double> ,double> > InteracaoB,InteracaoN,InteracaoS;
vector< vector< pair< pair<double,double> ,double> > > EscN;
double F(double d, double a) // função para encontrar as raizes
{
	return (a * exp(d)) - (4 * (d * d));
}
 
double dF(double d, double a) // derivada da função que queremos encontrar as raizes
{
	return a * exp(d) - 8 * d;	
}
 
void buildSigns(double a) // construir "tabela" de sinais
{
	bool sig = (F(-10, a) > 0);
	double last = -10.0;
	for(double i = -10; i <= 10; i += 0.5)
	{
		double Fx = F(i, a);
		sign[i + C] = (Fx > 0);
		if(sign[i + C] != sig)
		{
		    // adicionar no vector values os intervalos em que devemos procurar as raizes
			values.push_back({last, i}); // F(last) * F(i) < 0, respeitando o teorema do sinal
			sig = !sig;
		}
		last = i;
	}
}
 
 
// metodo da bisecao
double bisection(double l, double r, double a)
{
	InteracaoB.clear();
	contador = 0;
	double middle = -1, Fx = -1;
	bool sig = (F(l, a) > 0);
	while(fabs(r - l) > EPS)
	{
		middle = (l + r) / 2.0; 
		Fx = F(middle, a);
		if(!sig)
		{
			if(Fx < 0)		
				l = middle;
			else
				r = middle;
		}
		else
		{
			if(Fx > 0)
				l = middle;
			else
				r = middle;
		}
		InteracaoB.push_back({{middle,Fx},r-l});
		contador++;
		if(contador == iteracoes) break;
	}
	return middle;
}
 
// metodo de newton-raphson
double newtonRaphson(double xi, double a)
{
	contador = 0;
	double h = F(xi, a) / dF(xi, a);
	while(fabs(h) >= EPS)
	{
		h = F(xi, a) / dF(xi, a);
		xi -= h;
		InteracaoN.push_back({{xi,F(xi,a)},h});
		contador++;
		if(contador == iteracoes) break;
	}
	return xi;
}
 
// metodo da secante
double secant(double l,double r,double a)
{
	InteracaoS.clear();
	contador = 0;
	double fim, ini, c;
	do
	{
		ini = (l * F(r, a) - r * F(l, a)) / (F(r, a) - F(l, a));
		c = F(l, a) * F(ini, a);
		l = r;
		r = ini;
		if(fabs(c) < EPS)
			break;
		fim = (l * F(r, a) - r * F(l, a)) / (F(r, a) - F(l, a));
		InteracaoS.push_back({{r,F(r,a)},fim-ini});
	    contador++;
		if(contador == iteracoes) break;
	}
	while(fabs(fim - ini) >= EPS);
	return ini;
}
 
 
 
// as funções abaixo servem para imprimir os resultados
 
void printRootsWithBisection(double a)
{	
	puts("Roots - Bisection");
	for(auto p : values){
		double root = bisection(p.first, p.second, a);
		cout << fixed << setprecision(casas) << "\nraiz: " << root << '\n';
		puts("Interacoes:");
		puts("");
		cout << setw(casas/2+3) << "x" << setw(casas+casas/2+2) << "F(x)" << setw(0.75*casas+11) << "diferenca\n";
		for(auto i : InteracaoB)
			cout << fixed << setprecision(casas) << i.x << "  "<<((i.fx<0)?"":" ") << i.fx << "  "<<((i.dif<0)?"":" ") << i.dif << '\n';
	}
}
 
void printRootsWithNewtonRaphson(double a)
{
 
	puts("Roots - Newton Raphson");
	set<string> root;
	for(auto p : values)
	{
		string aux = to_string(newtonRaphson(p.first, a));
		if(root.find(aux)==root.end())
			EscN.push_back(InteracaoN);
		InteracaoN.clear();
		root.insert(aux);
		aux = to_string(newtonRaphson(p.second, a));
		if(root.find(aux)==root.end())
			EscN.push_back(InteracaoN);
		InteracaoN.clear();
		root.insert(aux);
	}
	int s=0;
	for(auto it : root){
		printf("\nraiz : %s\n", it.c_str());
		puts("Interacoes:");
		puts("");
		cout << setw(casas/2+3) << "x" << setw(casas+casas/2+2) << "F(x)" << setw(0.75*casas+11) << "diferenca\n";
		for(auto i:EscN[s])
			cout << fixed << setprecision(casas) << i.x << "  "<<((i.fx<0)?"":" ") << i.fx << "  "<<((i.dif<0)?"":" ") << i.dif << '\n';
		s++;
	}
}
 
void printRootsWithSecant(double a)
{	
	puts("Roots - Secant");
	for(auto it : values){
		double root = secant(it.first, it.second, a);
		cout << fixed << setprecision(casas) << "\nraiz: " << root << '\n';
		puts("Interacoes:");
		puts("");
		cout << setw(casas/2+3) << "x" << setw(casas+casas/2+2) << "F(x)" << setw(0.75*casas+11) << "diferenca\n";
		for(auto i : InteracaoS)
			cout << fixed << setprecision(casas) << i.x << "  "<<((i.fx<0)?"":" ") << i.fx << "  "<<((i.dif<0)?"":" ") << i.dif << '\n';
	}
}
 
int main()
{
	double a;
 
	puts("quantidade de iteracoes, precisao e o valor de \"a\" da funcao.");
	scanf(" %d %lf %lf", &iteracoes, &EPS, &a);
 
    double aux = EPS;
    while(int(aux) == 0)
        aux *= 10, casas++;
 	
	buildSigns(a);
 
	printRootsWithBisection(a);
	puts("\n");
 	printRootsWithNewtonRaphson(a);
 	puts("\n");
 	printRootsWithSecant(a);
 
	return 0;
}
