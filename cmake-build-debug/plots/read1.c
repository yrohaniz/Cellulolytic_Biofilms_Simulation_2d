/* This program reads the data from the cellulolytic_biofilms_simulation output
 and parameters from params to produce especially formatted files where the data 
 is organized in rows and columns. The formatted data is represented in the matrix
 form. The first number in the (1,1) place is the total number of data 
 points. The number of files is determined before creating them. This is acheived
 by reading the sim file and counting the number of time points. The first row in
 the formatted files represents the x-coordinate of the data points and the first 
 column represents the y-coordinate. This type of formatting is compatible with
 the gnuplot requirements to produce a heat-map-type plot.*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>

//****************************************************************FILE SCOPE INITIALIZATIONS*********************************************************|                                                        |
//**************************************************************************************************************************************************|
double DOM_LEN; // The length of the simulated domain                                                       
long NUM; // The grid refinement in both x and y direction. It is assumed that these are equal.                                 
const int NUM_CHARS=100; //Initialize the size of string variables used below.    
//**************************************************************************************************************************************************|
//**************************************************************************************************************************************************|

//This function goes through the white-spaces in a line of params file and grabs the value of the parameter needed for formatting
void param_iteration(long n, char *s, FILE *f) {for (int i=0; i<=n; i++) fscanf(f, "%s", s);}

int main() {
	int i, j, k, ntimes, count=0; //Declare ntimes for recording the number of time points
	long n; // The number of white-spaces in the params file
	char *ptr; // This pointer is used in strtol and is set the next character in the string stored by strtol
	FILE *f_in, *f_out;
	char buff[NUM_CHARS]; // Declaration of a string that stores the input data point
	char SIM_FILE[NUM_CHARS]; //Declaration of a string that stores the input data-file name
	double x, y; // x and y values of the data points
	double dx, dy; // increments along the x and y axis
	double *times, *tmp; // *times is used to store the time point values for further processing, *tmp is used to store the data corresponding to each time point
	char str[NUM_CHARS]; // Declaration of an string to arbitrary length
	char src[NUM_CHARS]; // ditto as above
	char sg[NUM_CHARS];// ditto as above
	struct stat st; // struct variable used for checking whether the folder where the files are saved, exists
	
	f_in = fopen ("params", "r");
	
	if (fgets(str, NUM_CHARS, f_in) !=NULL) { // Check whether params is opened and have contents
		fgets(str, NUM_CHARS, f_in); // Read the first line which is the number of white spaces
		n = strtol(str, &ptr, 10); //strtol reads the string in the first line and converts it to unsigned long integer
		fgets(str, NUM_CHARS, f_in); // Read the next line which is number of parameters
		printf("%s", str); // Print to the screen the number of parameters
		param_iteration(n,str,f_in); // go through white spaces within the line
		DOM_LEN = strtod(str, &ptr); // Pick the value of the domain length
		printf("domain length: %2.1f\n", DOM_LEN);
		param_iteration(n,str,f_in);
		NUM = strtol(str, &ptr, 10); // Pick the value of the grid refinement
		printf("refinement: %ld\n", NUM);
		param_iteration(n,str,f_in);
		strcpy(SIM_FILE, str); // Read the name of the input file
		printf("data file: %s\n", SIM_FILE);
		fscanf(f_in, "%s", str);
		if (feof(f_in)) { // Check whether the end of file is true and if not exit and print warning to screen
			printf("Parameter reading complete.\n");
		}
		else { 
			printf("ERROR in file!CHECK params!\n");
			exit(1);
		}
	}
	
	fclose(f_in); // Close the params file
	
	f_in = fopen (SIM_FILE, "r"); //Open the raw simulation file which contains all the times and corresponding biomss/carb values
	
	if ((times = (double *) malloc((NUM*NUM)*sizeof(double))) == NULL) {
        fprintf(stderr,"malloc failed\n");
        exit(1);
    }// Allocate memory to *times. Size is not known beforehand hence a big number is chosen
    if ((tmp = (double *) malloc((NUM*NUM)*sizeof(double))) == NULL) {
        fprintf(stderr,"malloc failed\n");
        exit(1);
    }// Allocate memory to *tmp. Size here is known which is equal to the total number of data points at each time point
    
    ntimes = 0; // Initialize ntimes to zero, its final value is determined by the while-loop
    
	while(1) {
		fscanf(f_in, "%s", buff); // Read the first data in input file (which is time). It is discarded by updating buff in the following loop
		if (feof(f_in)) {
			break;
		}// If end of file is reached break the loop
		times[ntimes] = strtod(buff, &ptr); // The first entry of a line in the sim file is the time point
		for (i=0; i<(NUM*NUM); i++) {
			fscanf(f_in, "%s", buff);
		} //The rest of the data in a particular line are read to the end of line
		ntimes++;
	}
	rewind(f_in); // Go to the beginning of the input_file
	times = (double *) realloc(times,ntimes*sizeof(double)); // Reallocate memory used to store time points to the true number of time points
	
	sprintf(src, "formatted_files_%ldx%ld/", NUM, NUM);
	
	if (stat(src, &st) == -1) {
    		mkdir(src, 0700);
	}// This if-statement checks the existence of the directory '/formatted_files' and creates it if not found
	
	// The following for-loop creates files for each time point and stores the corresponding sim data according to the desired format
	for (i=0; i<ntimes; i++) {
		dx = DOM_LEN/ (double) NUM; // Initialize the value of the increment for x values
		dy = DOM_LEN/ (double) NUM; // Initialize the value of the increment for y values
		x = dx/2.0; // Initialize the first x value to the abscissa of the center of cell
		y = dy/2.0; // Initialize the first y value to the ordinate of the center of cell
    	sprintf(str, "%s_formatted_t%d", SIM_FILE, (int) times[i]); // Append to the name of the file the corresponding value of the time point and store the result in str
    	//strcat(src, "/"); // Copy the string in quotes to src
    	strcpy(sg, src);
    	strcat(sg, str); // String concatenation i.e. append str to scr
    	f_out = fopen (sg, "w"); // After all the file name creation operations open/create the file
    	fprintf(f_out, "%ld ", NUM*NUM); // Save the first data to the output file which is the total number of points
		for (j=0; j<NUM; j++) {
			fprintf(f_out, "%f ", x); // Print the x-values in the first row
			x += dx; // Increment x
		} // This for-loop saves the x cooardinates to the file (after the total number of points from above) in the first row
		fprintf(f_out, "\n"); // Add newline after the loop is finished
		
		fscanf(f_in, "%s", buff); // Read the first data in input file (which is time). It is discarded by updating buff in the following loop
		for (j=0; j<NUM*NUM; j++) {
			fscanf(f_in, "%s", buff);
			tmp[j] = strtod(buff, &ptr);
		}// This for-loop stores the data corresponding to the current time point in *tmp
		
		for (j=0; j<NUM; j++) {
			fprintf(f_out, "%f ", y); // Print the y-value in the first column
			y += dy; // increment y
			for (k=0; k<NUM*NUM; k++) {
				if (k % NUM == j) {
					fprintf(f_out, "%f ", tmp[k]); // If the data belongs to the same y-value but different x-value save on the formatted file on the same line
				}// This for-loop goes through all the data in the tmp array and prints the ones that have the same y-value to the output file
			}
			fprintf(f_out, "\n"); // Add newline each time loop through j is performed (go to the next row)
		}
    	fclose(f_out); // Before creating the next formatted file close the current one
    }
	
	fclose(f_in); // Close the input file
	free(times);
	free(tmp);
	return 0;
}
