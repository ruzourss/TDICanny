#include <C_Matrix.hpp>
#include <C_Image.hpp>
#include <C_Arguments.hpp>
#include <C_File.hpp>
#include <math.h>
#include <iostream>
class ImageTDI
{
public:
	C_Arguments* argumentos;
	C_Image image;
	C_Image imageProcess;
	ImageTDI();
	~ImageTDI();

	C_Image imagenOriginal, imagenGauss, imagenGradiente,
		imagenNoMax, imagenBordeFuerte, imagenBordeDebil, 
		imagenOrientacion,imagenBinario, imagenFinal, contornos;
	C_Matrix kernelGauss;

	float t1,t2;

	/*
	Cargamos la imagen original
	*/
	void loadImage(const char* path) {
		imagenOriginal.ReadBMP(path);
		if (!imagenOriginal.Fail()) {
			imagenOriginal.WriteBMP("imagenesProcesadas/imagenOriginal.bmp");
		}	
	}


	/*
	comprobamos de que el fichero que estamos buscado existe
	*/
	bool checkFile(int arg, char **args) {
		argumentos = new C_Arguments(arg, args, 1, 2, "f", false);
		if (argumentos->Fail())
		{
			return false;
		}
		else {
			return (!C_FileExists(argumentos->Param(1))) ? false : true;
		}
	}

	const char* getPath() {
		return argumentos->Param(1);
	}

	void pause() {
		getchar();
	}
	/*
		Generamos la matriz gaussiana
	*/
	void generateKernelGaussian(int size, double sigma) {
		kernelGauss = C_Matrix(0,size-1,0,size-1);
		getGaussianKernel(sigma);
	}

	/*
	 Obtenemos un histograma sencillo de la imágen en un vector desde 0 a 255
	*/
	int getHistogram(int argc, char **argv) {
		C_Image image;
		//Comprobamos de que el fichero exista
		if (checkFile(argc, argv)) {
			image.ReadBMP(getPath());
			int histograma[255] = {};
			long row, col;
			for (row = image.FirstRow();row <= image.LastRow();row++) {
				for (col = image.FirstCol();col <= image.LastCol();col++) {
					histograma[(int)image(row, col)]++;
				}
			}
			int index;
			for (index = 0;index < 255;index++) {
				printf("%i=%i\n", index, histograma[index]);
			}
		}
		else {
			printf("File not found");
		}
		getchar();
		image.Free();
		return 0;
	}

	/*
		Calculamos el gradiente y guardamos el valor en una matriz nueva
	*/
	void applyGradient() {

		imagenGradiente = imagenOriginal;
		imagenGradiente.SetValue(0);
		imagenOrientacion = imagenGradiente;

		long width = imagenGradiente.LastCol();
		long height = imagenGradiente.LastRow();
		long startX = imagenGradiente.FirstCol();
		long startY = imagenGradiente.FirstRow();

		long deltaX, deltaY;
		double Gx, Gy, div;
		deltaX = deltaY = 1;
		long maxX = width - deltaX;
		long maxY = height - deltaY;
		double orientation = 0.0;
		for (long i = startY;i <= height;i++) {
			for (long j = startX;j <= width;j++) {
				Gx = Gy = 0.0;
				Gx = (((j < maxX) ? imagenOriginal(i, j + deltaX) : 0) - ((j >deltaX) ? imagenOriginal(i, j - deltaX) : 0)) / 2 * deltaX;
				Gy = (((i < maxY) ? imagenOriginal(i + deltaY, j) : 0) - ((i > deltaY) ? imagenOriginal(i - deltaY, j) : 0)) / 2 * deltaY;

				imagenGradiente(i, j) = sqrt(Gx * Gx + Gy * Gy);

				if (Gx == 0) {
					orientation = (Gy == 0) ? 0 : 90;
				}
				else
				{
					div = (double)abs((Gy / Gx));
					orientation = atan(div);

					if (orientation < 22.5) {
						orientation = 0;
					}
					else if (orientation < 67.5) {
						orientation = 45;
					}
					else if (orientation < 112.5) {
						orientation = 90;
					}
					else if (orientation < 157.5) {
						orientation = 135;
					}
					else {
						orientation = 0;
					}
				}

				imagenOrientacion(i, j) = orientation;
			}
		}

		imagenOrientacion.WriteBMP("imagenesProcesadas/orientacion.bmp");
		C_Image gradienteMasVisible(imagenGradiente);
		gradienteMasVisible.AddEscalar(50);
		gradienteMasVisible.WriteBMP("imagenesProcesadas/gradienteMagnitud.bmp");
	}
	/*
		Supresión no máxima
	*/
	void AplicarSupresionNoMaxima() {

		imagenNoMax = imagenGradiente;
		long width = imagenNoMax.LastCol();
		long height = imagenNoMax.LastRow();
		long startX = imagenNoMax.FirstCol();
		long startY = imagenNoMax.FirstRow();
		long margen = 5;
		double leftPixel = 0.0, rightPixel = 0.0;

		for (long i = startY;i <= height;i++) {
			for (long j = startX;j <= width;j++) {
				if ((i > startY + margen && i < height - margen) && (j > startX + margen && j < width - margen)) {
					switch ((int)imagenOrientacion(i, j))
					{
					case 0:
						leftPixel = imagenGradiente(i, j - 1);
						leftPixel = imagenGradiente(i, j + 1);
						break;
					case 45:
						leftPixel = imagenGradiente(i - 1, j + 1);
						leftPixel = imagenGradiente(i + 1, j - 1);
						break;
					case 90:
						leftPixel = imagenGradiente(i - 1, j);
						leftPixel = imagenGradiente(i + 1, j);
						break;
					case 135:
						leftPixel = imagenGradiente(i - 1, j - 1);
						leftPixel = imagenGradiente(i + 1, j + 1);
						break;
					}
					if (imagenGradiente(i, j) < leftPixel || imagenGradiente(i, j) < rightPixel) {
						imagenNoMax(i, j) = 0;
					}
				}
				else {
					imagenNoMax(i, j) = 0;
				}
			}
		}
		C_Image noMax(imagenNoMax);
		noMax.AddEscalar(40);
		noMax.WriteBMP("imagenesProcesadas/supresionNomaxima.bmp");
	}
	
	/*
		Aplicamos el filtro gaussiano para reducir el ruido
	*/
	void applyGaussian() {

		imagenGauss = imagenOriginal;

		long firstX = imagenGauss.FirstCol();
		long firstY = imagenGauss.FirstRow();

		long lastX = imagenGauss.LastCol();
		long lastY = imagenGauss.LastRow();

		long XLine, YLine;
		double gris, div;

		long firstXkernel = kernelGauss.FirstRow();
		long lastXKernel = kernelGauss.LastRow();

		int distance = (lastXKernel) / 2;

		for (long fila = firstY; fila <= lastY; fila++) {
			for (long columna = firstX; columna <= lastX; columna++) {
				gris = div = 0;
				for (int i = 0;i <= lastXKernel; i++) {
					XLine = fila + (i - distance);
					for (int j = 0;j <= lastXKernel; j++) {
						YLine = columna + (j - distance);
						if ((XLine >= firstY && XLine <= lastY) && (YLine >= firstX && YLine <= lastX)) {
							gris += kernelGauss(i, j) * imagenOriginal(XLine, YLine);
						}
						else {
							int r = fila + i - distance;
							int c = columna + j - distance;

							if (r < firstY) r = firstY;
							if (r > lastY) r = lastY;

							if (c < firstX) c = firstX;
							if (c > lastX) c = lastX;

							gris += kernelGauss(i, j) * imagenOriginal(r, c);
						}
					}
				}
				div = kernelGauss.Sum();
				if (div != 0) {
					gris /= div;
				}

				gris = (gris > 255) ? 255 : gris;
				gris = (gris < 0) ? 0 : gris;

				imagenGauss(fila, columna) = gris;
			}
		}
		imagenGauss.WriteBMP("imagenesProcesadas/Gaussiano.bmp");
	}
	
	void encontrarBordesFuertesYDebiles(int t1, int t2) {

		imagenBordeFuerte = imagenNoMax;
		imagenBordeFuerte.SetValue(0);
		imagenBordeDebil = imagenBordeFuerte;
		contornos = imagenBordeFuerte;

		int x = imagenOriginal.FirstRow();
		int y = imagenOriginal.FirstCol();

		int endX = imagenOriginal.LastRow();
		int endY = imagenOriginal.LastCol();

		
		double value = 0;

		for (int i = x;i < endX;i++) {
			for (int j = y; j < endY;j++) {
				value = imagenNoMax(i, j);
				if (value >= t2) {
					imagenBordeFuerte(i, j) = 255;
					contornos(i, j) = 1;
				}
				if ((value<t2) && value >= t1) {
					imagenBordeDebil(i,j) = 255;
					contornos(i, j) = 2;
				}
			}
		}
		
		imagenBordeFuerte.WriteBMP("imagenesProcesadas/bordesFuertes.bmp");
		imagenBordeDebil.WriteBMP("imagenesProcesadas/bordesDebiles.bmp");

	}

	void buscarEnlaces(int i, int j) {
		double valor;
		int k;
		int x[8] = { 1, 1, 0, -1, -1, -1, 0, 1 }, y[8] = { 0, 1, 1, 1, 0, -1, -1, -1 };
		for (k = 0; k<8; k++) {
			valor = contornos(i + x[k], j + y[k]);
			if (valor == 2) {
				contornos(i + x[k], j + y[k]) = 1;
				buscarEnlaces(i + x[k], j + y[k]);
			}
		}
	}

	void Histeresis() {
		int x = imagenNoMax.FirstRow();
		int y = imagenNoMax.FirstCol();
		int endX = imagenNoMax.LastRow();
		int endY = imagenNoMax.LastCol();
		
		for (int i = x; i<endX; i++) {
			for (int j = y; j<endY; j++) {
				if (i != x && j != y && i != endX && j != endY && contornos(i, j) == 1) {
					buscarEnlaces(i, j);
				}
			}
		}
		contornos.SetValue(2, 0);
		contornos.SetValue(1, 255);
		contornos.WriteBMP("imagenesProcesadas/imagenFinal.bmp");
	}

	void pintarBordes() {
		imagenFinal = imagenOriginal;
		imagenFinal.SetValue(0, 1);
		imagenFinal.palette(0, C_RED) = 0;
		imagenFinal.palette(0, C_GREEN) = 255;
		imagenFinal.palette(0, C_BLUE) = 0;
		
		int x = contornos.FirstRow();
		int y = contornos.FirstCol();
		int endX = contornos.LastRow();
		int endY = contornos.LastCol();
		for (int i = x; i<endX; i++) {
			for (int j = y; j<endY; j++) {
				if (contornos(i, j) == 255) {
					imagenFinal(i, j) = 0;
				}
			}
		}
		imagenFinal.WriteBMP("imagenesProcesadas/imagenFinalConBorde.bmp");
	}

	void freeMemory() {
		imagenOriginal.Free();
		imagenGauss.Free();
		imagenGradiente.Free();
		imagenNoMax.Free();
		imagenBordeFuerte.Free();
		imagenBordeDebil.Free();
		imagenOrientacion.Free();
		imagenFinal.Free();
		contornos.Free();
		kernelGauss.Free();
	}

private:
	double getGaussian2D(double sigma, int x, int y) {
		return exp(-((x * x + y * y) / (2 * sigma * sigma)));
	}

	void getGaussianKernel(double sigma) {

		int lastRow = kernelGauss.LastRow();

		int range = lastRow / 2;
		double sum = 0;

		for (int i = 0, y = -range;i <= lastRow;i++, y++) {
			for (int j = 0, x = -range;j <= lastRow;j++, x++) {
				kernelGauss(i, j) = getGaussian2D(sigma, x, y);
			}
		}
		sum = kernelGauss.Sum();
		kernelGauss.DivideEscalar(sum);
	}
};

ImageTDI::ImageTDI()
{
}

ImageTDI::~ImageTDI()
{
}
int main(int argc, char **argv)
{
	ImageTDI* imagen = new ImageTDI();
		
		C_Print("*********************************************************");
		C_Print("*\tDETECCION DE BORDES CON EL METODO DE CANNY\t*");
		C_Print("*\tTautvydas Bagocius\t\t\t\t*");
		C_Print("*********************************************************");

		char nombreImagen[100];
		int sizeMascara, sigma, t1, t2;

		do {
			imagen->imagenOriginal.Clear();
			cout << "Introduce el nombre la imagen: ";
			cin.getline(nombreImagen,100);
			imagen->loadImage(nombreImagen);
		} while (imagen->imagenOriginal.Fail());
		C_Print("Cargando imagen....");

		do {
			cout << "Introduce el tamaño de la mascara: ";
			cin >> sizeMascara;
			if (sizeMascara <= 0 || sizeMascara > 10) {
				printf("Rango de la mascara debe estar entre 1-10\n");
				system("pause");
				exit(0);
			}
		} while (sizeMascara <= 0 || sizeMascara > 10);

		do {
			cout << "Introduce la desviacion estadar: ";
			cin >> sigma;
			if (sigma <= 0 || sigma > 10) {
				printf("Rango de la desviacion estadar debe estar entre 1-10\n");
				system("pause");
				exit(0);
			}
		} while (sigma <= 0 || sigma > 10);

		do {
			cout << "Introduce el umbral minimo: ";
			cin >> t1;
			if (t1 < 0 || t1 > 255) {
				printf("Rango del umbral mínimo debe estar entre 0-255\n");
				system("pause");
				exit(0);
			}
		} while (t1 < 0 || t1 > 255);

		do {
			cout << "Introduce el umbral maximo: ";
			cin >> t2;
			if (t2 < 0 || t2 > 255) {
				printf("Rango del umbral maximo debe estar entre 0-255\n");
				system("pause");
				exit(0);
			}
		} while (t2 < 0 || t2 > 255);
		 C_Print("Aplicando filtros.... ");
		 C_Print("*********************************************************");
		 imagen->generateKernelGaussian(sizeMascara,sigma);
		 imagen->applyGaussian();
		 printf("Termino el filtro Gaussiano\n");
		 imagen->applyGradient();
		 printf("Gradiente aplicado\n");
		 imagen->AplicarSupresionNoMaxima();
		 printf("Supresion no maxima aplicado\n");
		 imagen->encontrarBordesFuertesYDebiles(t1,t2);
		 printf("Buscanado bordes\n");
		 imagen->Histeresis();
		 printf("Histeresis aplicado\n");
		 imagen->pintarBordes();
		 printf("Bordes pintados\n");
		 C_Print("*********************************************************");
		 imagen->freeMemory();
		 system("pause");
}

