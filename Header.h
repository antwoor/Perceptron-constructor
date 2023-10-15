#pragma once
#include <vector>
#include <ctime>
#include <math.h>
using namespace std;
class Way
{
public: 
	Way() {};
	vector<double> v;  
	int l;   
	Way(double arr[], int length)
	{
		l = length;
		v.resize(l);
		for (int i = 0; i < l; i++)
		{
			v[i] = arr[i];
		}
	}
	Way(int l)      
	{
		this->l = l;    
		v.resize(l);     
	}
	Way(vector<double> values)
	{
		l = values.size();
		v.resize(l);
		for (int i = 0; i < l; i++)
		{
			v[i] = values[i];
		}
	}
	void ShowWay() {
		for (int i = 0; i < l; i++)
		{
			cout << v[i] <<" ";
		}
	}
	double GetWay(int i)
	{
		return v[i];
	}
	double SetWay(int i, double value)
	{
		v[i] = value;
		return 0;
	}
	~Way() {};
};

class matrix
{
public:
	double** a;
	int m, n;
	matrix() {};
	matrix(int m, int n)
	{
		this->m = m;        
		this->n = n;
		a = new double* [m];          
		for (int i = 0; i < m; i++)
		{
			a[i] = new double[n];                                       
			for (int j = 0; j < n; j++)
			{
				a[i][j] = 0;
			}
		}
	}
	matrix(int m, int n, double RandNum)
	{
		srand(time(NULL));
		this->m = m;        
		this->n = n;
		a = new double* [m];          
		for (int i = 0; i < m; i++)
		{
			a[i] = new double[n];                                       
			for (int j = 0; j < n; j++)
			{
				RandNum =rand() %10;
				RandNum = RandNum / 10 + 1;
				a[i][j] = RandNum;
			}
		}
	}
public:
	double GetMatrix(int i, int j)
	{
		return a[i][j];
	}
public:

	double SetMatrix(int i, int j, double value)
	{
		return a[i][j] = value;
	}
	void ShowMatrix()
	{
		for (int i = 0; i < m; i++)
		{
			for (int j = 0; j < n; j++)
			{
				cout << a[i][j] << " ";
			}
			cout <<  endl;
		}
	}
	matrix transpond(matrix a)
	{
		int atn, atm;
		atn = a.m;
		atm = a.n;
		matrix at(atm, atn);
		for (int i = 0; i < atm; i++)
		{
			for (int j = 0; j < atn; j++)
			{
				at.SetMatrix(i, j, a.GetMatrix(j, i));
			}
		}
		return at;
	}
	~matrix() {};
};

class Network
{
public:
	Network() {};
	struct GLayer
	{
		public:
			Way x, z, df;            
	};
	Way* deltas = new Way[2];     
	GLayer* cur = new GLayer[2];     
	vector<matrix> weights;      
	vector<matrix> DeltasW;    
	int layersN;
	Network(vector<int> sizes)         
	{
		srand(time(NULL));
		double random = rand();
		layersN = sizes.size() - 1;     
		weights.resize(layersN);         
		DeltasW.resize(layersN);     
		cur = new GLayer[layersN];                  
		deltas = new Way[layersN];     
		for (int i = 0; i < layersN; i++)
		{
			weights[i] = matrix(sizes[i], sizes[i+1], random);   
			DeltasW[i] = matrix(sizes[i], sizes[i + 1]);
			cur[i].x = Way(sizes[i]);    
			cur[i].z = Way(sizes[i+1]);     
			cur[i].df = Way(sizes[i+1]);     
			deltas[i] = Way(sizes[i+1]);     
		}
	}	                 
		  Way Forward(Way in)
		  {
			  for (int i = 0; i < layersN; i++)
			  {
				  if (i == 0)
				  {
					  for (int j = 0; j < in.l; j++)
					  {
						  cur[i].x.SetWay(j, in.GetWay(j));
					  }
				  }
				  else
				  {
					  for (int j = 0; j < cur[i-1].z.l; j++)
					  {
						  cur[i].x.SetWay(j, cur[i - 1].z.GetWay(j));
					  }
				  }
				  for (int k = 0; k < weights[i].n; k++)  
				  {
					  double y = 0; 
					  for (int j = 0; j < weights[i].m; j++)
					  {
						  y += weights[i].GetMatrix(j, k) * cur[i].x.GetWay(j);
						  cur[i].z.SetWay(k, (1 / (1 + exp(-y)))); 
						  cur[i].df.SetWay(k, ((cur[i].z.GetWay(k)) * (1 - cur[i].z.GetWay(k))));        
					  }
				  }
			  }
			  return cur[layersN-1].z;  
		  }
public:
	void BackProp(Way output, double& error)
	{    
		int last = layersN-1;      
		error = 0;   
		for (int i = 0; i < output.l; i++)
		{
			double dif = output.GetWay(i) - cur[last].z.GetWay(i);        
			deltas[last].SetWay(i, (dif * cur[last].df.GetWay(i)));    
			error += (dif * dif);        
		}
		error = error / output.l;        
		for (int k = last; k > 0 ; k--)   
		{
			for (int i = 0; i < weights[k].m; i++) 
			{
				deltas[k - 1].SetWay(i, 0);      
				for (int j = 0; j < weights[k].n; j++)
				{
					deltas[k - 1].SetWay(i, (deltas[k - 1].GetWay(i) + weights[k].GetMatrix(i, j) * deltas[k].GetWay(j)));
					deltas[k - 1].SetWay(i, deltas[k - 1].GetWay(i) * cur[k - 1].df.GetWay(i));
				}
			}
		}
	}
	void UpdWeights(double alpha, double E)
	{
		for (int k = 0; k < layersN; k++)
		{
			for (int i = 0; i < weights[k].m; i++)
			{
				for (int j = 0; j < weights[k].n; j++)
				{
					double gradient = deltas[k].GetWay(j) * cur[k].x.GetWay(i);
					double deltaweight = E * gradient  +alpha * DeltasW[k].GetMatrix(i, j);
					DeltasW[k].SetMatrix(i, j, deltaweight);
					weights[k].SetMatrix(i, j, (weights[k].GetMatrix(i, j) + deltaweight));                   
				}
			}
		}
	}
	public: 
		void ShowNetwork()
		{
			for (int i = 0; i < layersN; i++)
			{
				cur[i].x.ShowWay();
				cout << endl;
				cout << "||||||||||||" << endl;
				weights[i].ShowMatrix();
				cout << "////////////" << endl;
			}
		};
		void Train(vector<Way> X, vector<Way> Y, double alpha, double eps, int epochs, double E)
		{
			int epoch = 1;   
			double error = 0;  
			do
			{
				for (int i = 0; i < X.size(); i++)
				{
					Forward(X[i]);
					BackProp(Y[i], error);
					UpdWeights(alpha, E);
				}
				cout << "epoch: ";
				cout <<epoch;
				cout << " error: ";
				cout<< error;
				cout << " " << endl;
				epoch++;
			} while (epoch <epochs+1);
		}
		~Network() {};
};

