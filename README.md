# CACS_Project
Project for the course of Cryptography and Architectures for Computer Security, by Davide Li Calsi and Giorgio Romeo.

## Usage

### Loading a generator and parity-check matrix
Copy the txt file that you can find at the url https://decodingchallenge.org/Challenges/SD/SD_xxx_0 ,
where xxx represents the parameter n, and paste it in the file Utilities/info.txt.

### Finding the best parameters
To compute the best parameters for decoding, i.e. the parameters that minimize complexity, run

```sh
python find_param.py n
```

where n must be equal to the size that you set in the info.txt file. The script assumes a code-rate of R=0.5, since this is the default one used in the challenge website. Then, simply copy them in the main.c file.

### Running
Compile with a simple

```sh
gcc -o main main.c -lm -ggdb3
```

or using 

```sh
make
```

Then run with

```sh
./main
```

The code will automatically generate a code and an error pattern and attempt to decode it. It will print the
guessed codewords and their Hamming Distance from the vector to decode as they are generated.