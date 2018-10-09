#include <iostream>
using namespace std;

int main() {
	int a[10], swap;
	for (int i = 0; i < 10; ++i)
		cin >> a[i];
	for (int i = 0; i < 10; ++i) {
		for (int j = i; j < 10; ++j){
			if (a[j] > a[i]){
				swap = a[i];
				a[i] = a[j];
				a[j] = swap;
			}
		}
	}

	//for (int i = 0; i < 10; ++i)
		//cout << a[i] << " ";
	//return 0;
}