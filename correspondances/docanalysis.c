// I N C L U D E S ///////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <string.h>
#include <ctype.h>

// D E F I N E S /////////////////////////////////////////////////////////////
 
#define MAX_WORD_LENGTH		64
#define MAX_FILE_NAME       256

// T Y P E S /////////////////////////////////////////////////////////////////

typedef enum { FALSE = 0, TRUE = 1 } boolean_t;

typedef unsigned long dword_t;

// S T R U C T U R E S ///////////////////////////////////////////////////////

typedef struct node_s
{
	char	word[MAX_WORD_LENGTH];
	dword_t	frequency;
	struct	node_s* next;
} node_t;

typedef struct
{
	node_t*			head;
	node_t*			tail;
	node_t			**words_ptrs;
	dword_t			size;
} word_list_t;

// G L O B A L S //////////////////////////////////////////////////////////////


// F U N C T I O N S //////////////////////////////////////////////////////////

void InitWordList(word_list_t *word_list)
{
	word_list->head		  = NULL;
	word_list->tail		  = NULL;
	word_list->size       = 0;
	word_list->words_ptrs = NULL;
} // end WordListInit

///////////////////////////////////////////////////////////////////////////////

void AddToWordList(word_list_t *word_list, const char *word)
{
	node_t *new_node;

	new_node = (node_t*)malloc(sizeof(node_t));

	strcpy(new_node->word, word);

	new_node->frequency = 1;
	new_node->next      = NULL;
	
	if (word_list->head == NULL)
	{
		word_list->head = new_node;
		word_list->tail = word_list->head;
	} // end if
	else
	{
		word_list->tail->next = new_node;
		word_list->tail = word_list->tail->next;
	} // end 

	word_list->size++;
} // end AddToWordList


///////////////////////////////////////////////////////////////////////////////

boolean_t InWordList(word_list_t *word_list, char* word)
{
	node_t *p;

	if (word_list->head == NULL )
		return FALSE;

	for (p = word_list->head; p != NULL; p = p->next)
	{
		if (strcmp(p->word, word) == 0)
		{
			p->frequency++;
			return TRUE;
		}
	}

	return FALSE;
} // end InWordList

///////////////////////////////////////////////////////////////////////////////

void DestroyWordList(word_list_t *word_list)
{
	node_t *p, *prev;

	p = word_list->head;
	while (p != NULL)
	{
		prev = p;
		p = p->next;
		if (prev)
		{
			free(prev);
			prev = NULL;
		}
	} // end while
	
	if (word_list->words_ptrs)
	{
		free(word_list->words_ptrs);
		word_list->words_ptrs = NULL;
	}

	InitWordList(word_list);

} // end DestroyWordList

///////////////////////////////////////////////////////////////////////////////

int CompareWords( const void *arg1, const void *arg2 )
{
	node_t  **p1, **p2, *node_1, *node_2;

	p1 = (node_t**)arg1;
	p2 = (node_t**)arg2;

	node_1 = (node_t*)(*p1);
	node_2 = (node_t*)(*p2);

		
	if (node_1->frequency < node_2->frequency)
		return (1);
	else
	if (node_1->frequency > node_2->frequency)
		return (-1);
	else
		return strcmp( node_1->word, node_2->word );

} // end CompareWords

///////////////////////////////////////////////////////////////////////////////

void SortWordList(word_list_t *word_list)
{
	node_t *p;
	int count;

	if (!word_list->head)
		return;

	word_list->words_ptrs = (node_t**)malloc(sizeof(node_t*)*word_list->size);

	count = 0;
	for (p = word_list->head; p != NULL; p = p->next)
	{
		word_list->words_ptrs[count++]	= p;
	}
	
	qsort((void*)word_list->words_ptrs, (size_t)word_list->size, sizeof(node_t*), CompareWords );
} // end SortWordList

///////////////////////////////////////////////////////////////////////////////

void SaveSortedWordListToFile(FILE *fp, word_list_t *word_list)
{
	dword_t i;

	if (!word_list->words_ptrs)
		return;

	for (i=0; i<word_list->size; i++)
	{
		fprintf(fp, "%s %d\n", word_list->words_ptrs[i]->word,
			                   word_list->words_ptrs[i]->frequency);
	}
} // end SaveSortedWordListToFile

///////////////////////////////////////////////////////////////////////////////

void StrToLower(char* str)
{
	char *p = NULL;

	for( p = str; p < str + strlen(str); p++ )
	{
		*p = tolower(*p);
	}
	
} // end StrToLower

///////////////////////////////////////////////////////////////////////////////

boolean_t InsertWordsToList(word_list_t *word_list, char *buffer, dword_t num_chars)
{

	char ch, word[MAX_WORD_LENGTH];

	boolean_t word_found;

	dword_t i, word_num_chars;


	if (!buffer)
	{
		return FALSE;
	}

	word_found = FALSE;

	for (i=0; i<num_chars; i++)
	{
		ch = buffer[i];

		if ( !word_found )
		{
			if (isalpha(ch))
			{
				word_found = TRUE;
				word_num_chars = 0;
				word[word_num_chars++] = ch;
			} // end 
		} // end if
		else
		{
			if (!isalpha(ch))// && (ch != '-'))
			{
				word_found = FALSE;
				if ( word_num_chars > 1 )
				{						
					word[word_num_chars++] = '\0';

					StrToLower(word);

					if (!InWordList(word_list, word))
						AddToWordList(word_list, word);

				} // end if
			} // end if
			else
			{
				word[word_num_chars++] = ch;
			} // end else
		} // end else
	} // end for i


	return TRUE;
} // end InsertWordsToList


///////////////////////////////////////////////////////////////////////////////

int main(int argc, char* argv[])
{
	word_list_t word_list;

	char *buffer = NULL;

	dword_t file_size;

	FILE *fp;

	char file_in[MAX_FILE_NAME], 
		 file_out[MAX_FILE_NAME];


    if (argv[1])
    {
        strcpy(file_in, argv[1]);
    }
    else
    {
		printf("parameters: file_in [file_out]\n");
        return 0;
    }

    if (argv[2])
		strcpy(file_out, argv[2]);
    else
		strcpy(file_out, "docwords.txt");

	
	InitWordList(&word_list);

	if (!(fp = fopen(file_in, "rb")) )
    {
        printf("Error: cannot open input file!!!!\n");
		return 0;
    }

	// find file size
	fseek(fp, 0, SEEK_END);
	file_size = ftell(fp);
	fseek(fp, 0, SEEK_SET);


	// allocate memory buffer to store characters
	buffer = (char*)malloc(file_size);

	fread(buffer, sizeof(char), file_size, fp);
	fclose(fp);

	InsertWordsToList(&word_list, buffer, file_size);

	free(buffer);

	// Sort word list
	SortWordList(&word_list);


	if (!(fp = fopen(file_out, "w")) )
    {
        printf("Error: cannot open output file!!!!\n");
		return 0;
    }

	SaveSortedWordListToFile(fp, &word_list);

	fclose(fp);

	DestroyWordList(&word_list);

	printf("%s is created\n", file_out);

	return 0;
}  // end main
 
