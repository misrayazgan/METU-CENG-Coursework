#include <openssl/evp.h>
#include <iostream>
#include <fstream>

using namespace std;


/*********read the contents of a file into an unsigned char array************/
unsigned char* readFile(const string &fileName)
{
	FILE* file;
	unsigned char* content;
	long file_size;

	file = fopen(fileName.c_str(), "rb");

	fseek(file, 0, SEEK_END);		// go to the end of the file
	file_size = ftell(file);		// returns current position in stream (file size)
	
	rewind(file);					// set position of stream to the beginning

	content = (unsigned char*)malloc(sizeof(unsigned char) * file_size);
	fread(content, 1, file_size, file);

	fclose(file);
	return content;
}

/********compare two unsigned char arrays, return whether or not they are the same**********/
int compareArrays(unsigned char* a, unsigned char* b, int size)
{
	for(int i = 0; i < size; i++)
	{
		if(a[i] != b[i])
			return 0;
	}

	return 1;
}

/***********find how many words(lines) there are in the words.txt***********/
int findWordNumber(const string &fileName)
{
	ifstream file;
	int numberOfLines = 0;

	char* line;
	line = (char*)malloc(sizeof(char) * 30);

	file.open(fileName.c_str(), ifstream::in);
	
	while(!file.eof())
	{
		file.getline(line, 30);
		numberOfLines++;
	}

	file.close();
	free(line);

	return numberOfLines;
}

/**********read the given words into an array**********/
unsigned char** readWords(const string &fileName)
{
	FILE* file;
	file = fopen(fileName.c_str(), "rb");

	long file_size;

	unsigned char* content;

	fseek(file, 0, SEEK_END);
	file_size = ftell(file);
	rewind(file);

	content = (unsigned char*)malloc(sizeof(unsigned char) * file_size);

	fread(content, sizeof(unsigned char), file_size, file);

	fclose(file);

	int numberOfLines = findWordNumber(fileName);

	unsigned char** lines;
	lines = (unsigned char**)malloc(sizeof(unsigned char*) * numberOfLines);
	for(int i = 0; i < numberOfLines; i++)
	{
		*(lines + i) = (unsigned char*)malloc(sizeof(unsigned char) * 30);
	}

	/**********fill the lines array with ' '**********/
	for(int i = 0; i < numberOfLines; i++)
	{
		for(int j = 0; j < 30; j++)
		{
			lines[i][j] = ' ';
		}
	}

	/**********get words to the lines array***********/
	int k = 0;
	for(int i = 0; i < numberOfLines; i++)
	{
		for(int j = 0; j < 30; j++)
		{
			if(content[k] != '\n')
			{
				lines[i][j] = content[k++];
				
			}
			else
			{
				k++;
				break;
			}
		}
	}

	return lines;
}

/********AES-128-CBC encryption**********/
int encrypt(unsigned char* plaintextArray, unsigned char* givenCipher, unsigned char* key)
{
	int len1 = 0;
	int len2 = 0;

	int plaintextLength = 21;
	
	unsigned char* cipherArray;
	cipherArray = (unsigned char*)malloc(32 * sizeof(unsigned char));

	unsigned char iv[] = {'\0','\0','\0','\0','\0','\0','\0','\0','\0','\0','\0','\0','\0','\0','\0','\0'};

	unsigned char* mykey;
	mykey = (unsigned char*)malloc(16 * sizeof(unsigned char));

	int i = 0;
	while(key[i] != '\0')
	{
		mykey[i] = key[i];
		i++;
	}
	
	for(i; i < 16; i++)
	{
		mykey[i] = '\0';
	}

	// Create a new cipher context
	EVP_CIPHER_CTX *context;
	context = EVP_CIPHER_CTX_new();

	EVP_EncryptInit_ex(context, EVP_aes_128_cbc(), NULL, mykey, iv);	
	
	if(!EVP_EncryptUpdate(context, cipherArray, &len1, plaintextArray, plaintextLength))
	{
		EVP_CIPHER_CTX_free(context);
		return 0;
	}

	if(!EVP_EncryptFinal_ex(context, cipherArray + len1, &len2))
	{
		EVP_CIPHER_CTX_free(context);
		return 0;
	}

	int same = compareArrays(givenCipher, cipherArray, len1 + len2);
	
	EVP_CIPHER_CTX_free(context);

	return same;
}


int main()
{
	unsigned char* plaintextArray = readFile("plaintext.txt");

	unsigned char* givenCipher = readFile("ciphertext");

	unsigned char** words = (unsigned char**)readWords("words.txt");
	int numberOfWords = findWordNumber("words.txt");

	for(int i = 0; i < numberOfWords; i++)
	{
		unsigned char key[16];

		for(int j = 0; j < 16; j++)
		{
			if(words[i][j] == ' ' || j == 15)
			{
				key[j] = '\0';
				break;
			}
			else
			{
				key[j] = words[i][j];
			}
		}

		int same = encrypt(plaintextArray, givenCipher, key);

		if(same)
			cout << key << endl;
		
	}

}