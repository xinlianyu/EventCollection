#include<iostream.h>
#include<conio.h>

void main()
{
	clrscr();
	int a[10], n, i, j;
	int percent;
	int count;
	cout << "
		enter the size of array";
		cin >> n;
	cout << "
		enter the values";
		for (i = 0; i <= n - 1; i++)
		{
			cout << " a[" << i << "]" << "=";
			cin >> a[i];
		}

	for (i = 0; i <= n - 1; i++)
	{
		count = 0;
		for (j = 0; j <= n - 1; j++)
		{
			if (a[i]>a[j])
			{
				count = count + 1;
			}
		}
		percent = (count * 100) / (n - 1);
		cout << "

			the percentile of"<<"a["<<i<<"]"<<percent;
	}
	getch();

}