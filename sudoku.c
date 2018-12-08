#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <limits.h>
#include <string.h>
#include <stdbool.h>
#include <stddef.h>
#include <omp.h>
#include <pthread.h>
#include <mpi.h>

#define MAX_BDIM 8

enum SOLVE_STRATEGY {SUDOKU_SOLVE, SUDOKU_COUNT_SOLS};
#ifndef SUDOKU_SOLVE_STRATEGY
#define SUDOKU_SOLVE_STRATEGY SUDOKU_SOLVE
#endif

typedef struct sudoku
{
    unsigned char bdim;
    unsigned char dim;
    unsigned char peers_size;

    unsigned short solved;
    
    unsigned short **unit_list;
    unsigned short  **peers;
    unsigned long long *values;
    
    unsigned long long sol_count;

    int status;
} Sudoku;

static inline void cell_v_set(unsigned long long *v, int p);
static inline void cell_v_unset(unsigned long long *v, int p);
static inline int cell_v_get(unsigned long long *v, int p);
static inline int cell_v_count(unsigned long long *v);
static inline int digit_get (unsigned long long *v);

Sudoku* create_sudoku(unsigned char bdim, unsigned char *grid);
void init(Sudoku *s);
int parse_grid(Sudoku *s, unsigned char *grid);
int assign (Sudoku *s, int i, int d);
int eliminate (Sudoku *s, int i, int d);
void destroy_sudoku(Sudoku *s);
void display(Sudoku *s);
int task_search (Sudoku *s, int status);
Sudoku* copy_sudoku(Sudoku * s);
int distribution(Sudoku *s, int id, int ntasks);

bool jaDividiuProcessos = false, terminouLocalmente = false;
int world_size, world_rank, alguemTerminou = 0, flag[] = {0, 0, 0, 0, 0, 0};
pthread_t id;
MPI_Request polling;

void* initPolling (void* args) {

    MPI_Wait(&polling, MPI_STATUS_IGNORE);
    alguemTerminou = 1;
    return NULL;
}

static int search (Sudoku *s, int argMinI, int argK);

static inline void cell_v_set(unsigned long long *v, int p)
{
    *v |= (unsigned long long)1 << (p - 1);
}

static inline void cell_v_unset(unsigned long long *v, int p)
{
    *v &= ~((unsigned long long)1 << (p - 1));
}

static inline int cell_v_get(unsigned long long *v, int p)
{
    return !!(*v & ((unsigned long long)1 << (p - 1)));
}

static inline int cell_v_count(unsigned long long *v) 
{
    return __builtin_popcountll(*v);
}

static inline int digit_get (unsigned long long *v) 
{
    if (cell_v_count(v) != 1) return -1;
    return __builtin_ctzll(*v) + 1;
}

Sudoku* create_sudoku(unsigned char bdim, unsigned char *grid)
{
    assert(bdim <= MAX_BDIM);

    Sudoku *s = malloc(sizeof(Sudoku));
    
    s->bdim = bdim;
    unsigned char dim = bdim * bdim;
    s->dim = dim;
    s->peers_size = 3 * dim - 2 * bdim -1;
    s->sol_count = 0;
    s->solved = dim * dim;

    int i;
    s->unit_list = malloc(dim * dim * sizeof(unsigned short*));
    assert(s->unit_list);
    for( i = 0 ; i < dim * dim ; i++ )
    {
        s->unit_list[i] = calloc( 3 * dim, sizeof(unsigned short));
        assert(s->unit_list[i]);
    }

    s->peers = malloc(dim * dim * sizeof(unsigned short*));
    assert(s->peers);
    for( i = 0 ; i < dim * dim ; i++)
    {
        s->peers[i] = calloc(s->peers_size, sizeof(unsigned short));
        assert(s->peers[i]);
    }

    s->values = calloc(dim * dim, sizeof(unsigned long long));
    assert(s->values);
    
    init(s);
    if (!parse_grid(s, grid)) 
    {
        printf("Error parsing grid\n");
        destroy_sudoku(s);
        return 0;
    }
    return s;
}

void init(Sudoku *s)
{
    int i, j, l, k, pos;

    for( i = 0 ; i < s->dim ; i++ )
    {
        int ibase = i / s->bdim * s->bdim;

        for( j = 0 ; j < s->dim ; j++ )
        {
            for( pos = 0 ; pos < s->dim ; pos++ )
            {
                s->unit_list[i*s->dim + j][pos] = i*s->dim + pos;
                s->unit_list[i*s->dim + j][pos + s->dim] = pos*s->dim + j;
            }
            
            int jbase = j / s->bdim * s->bdim;
            for ( pos = 0, k = 0 ; k < s->bdim ; k++ )
            {
                for (l = 0 ; l < s->bdim ; l++ , pos++ )
                {
                    s->unit_list[i*s->dim + j][pos + 2*s->dim] = (ibase + k)*s->dim + jbase + l;
                }
            }
        }
    }

    for( i = 0 ; i < s->dim ; i++)
    {
        for( j = 0 ; j < s->dim ; j++)
        {
            pos = 0;

            for (k = 0; k < s->dim; k++) 
            {
                if (s->unit_list[i*s->dim + j][k] % s->dim != j)
                {
                    s->peers[i*s->dim + j][pos++] = s->unit_list[i*s->dim + j][k]; 
                }
            }

            for (k = 0; k < s->dim; k++) 
            { 
                int sq = s->unit_list[i*s->dim + j][k + s->dim];
                if (sq / s->dim != i)
                {
                    s->peers[i*s->dim + j][pos++] = sq;
                }

                sq = s->unit_list[i*s->dim + j][k + 2*s->dim];
                if (sq / s->dim != i && sq % s->dim != j)
                {
                    s->peers[i*s->dim + j][pos++] = sq;
                }
            }
        }
    }

    assert(pos == s->peers_size);
}

int parse_grid(Sudoku *s, unsigned char *grid)
{
    int i, k, dimxdim = s->dim * s->dim;
    for ( i = 0 ; i < dimxdim ; i++ )
    {
        for (k = 1; k <= s->dim; k++)
        {
            cell_v_set(&s->values[i], k);
        }
    }
    
    for ( i = 0 ; i < dimxdim ; i++)
    {
        if (grid[i] > 0 && !assign(s, i, grid[i]))
        {
            return 0;
        }
    }

    return 1;
}

int assign (Sudoku *s, int i, int d) 
{
    for (int d2 = 1; d2 <= s->dim; d2++)
    {
        if (d2 != d)
        {
            if (!eliminate(s, i, d2))
            {
               return 0;
            }
        }
    }

    return 1;
}

int eliminate (Sudoku *s, int i, int d) 
{
    int k, ii, cont, pos;
    
    if (!cell_v_get(&s->values[i], d))
    { 
        return 1;
    }

    cell_v_unset(&s->values[i], d);

    int count = cell_v_count(&s->values[i]);

    if (count == 0)
    {
        return 0;
    } 
    else if (count == 1)
    {
        s->solved--;
        for ( k = 0 ; k < s->peers_size ; k++ )
        {
            if (!eliminate(s, s->peers[i][k], digit_get(&s->values[i])))
            {
                return 0;
            }
        }
    }

    for ( k = 0 ; k < 3 * s->dim; k += s->dim )
    {
        cont = 0;
        pos = 0;
        for (ii = 0; ii < s->dim; ii++)
        {
            if (cell_v_get(&s->values[s->unit_list[i][k + ii]], d))
            {
                cont++;
                pos = k + ii;
            }
        }

        if (cont == 0)
        {
            return 0;
        }
        else if (cont == 1) 
        {
            if (!assign(s, s->unit_list[i][pos], d))
            {
                return 0;
            }
        }
    }

    return 1;
}

void destroy_sudoku(Sudoku *s)
{
    int i;
    for( i = 0 ; i < s->dim * s->dim ; i++)
    {
        free(s->unit_list[i]);
        free(s->peers[i]);
    }

    free(s->unit_list);
    free(s->peers);
    free(s->values);
    free(s);
}

void display(Sudoku *s)
{
    printf("\n\nsudoku - %d\n", s->bdim);
    for (int i = 0; i < s->dim; i++)
    {
        for (int j = 0; j < s->dim; j++)
        {
            printf("%d ", digit_get(&s->values[i*s->dim + j]));
        }
        printf("\n");
    }
}

int apresentarResultados (Sudoku * copia) {
    display(copia);
    MPI_Cancel(&polling);
    terminouLocalmente = true;
    alguemTerminou = 1;
    MPI_Send(&alguemTerminou, 1, MPI_INT, world_rank == 0 ? 1 : 0, 123, MPI_COMM_WORLD);
    return 1;
}

void encontrarSquareMenosPossibilidades (Sudoku *s, int *min, int *minI) {
    for (int i = 0; i < s->dim * s->dim; i++)  {
        int used = cell_v_count(&s->values[i]);
        if (used > 1 && used < *min) {
            *min = used;
            *minI = i;
        }
    }
}

void encontrarSquareNumProcessadores (Sudoku *s, int *min, int *minI) {
    int wanted = omp_get_num_procs(); 
    for (int i = 0; i < s->dim * s->dim; i++) {
        int used = cell_v_count(&s->values[i]);
        if (wanted == used) {
            *min = used;
            *minI = i;
            break;
        }
    }
}

Sudoku* copy_sudoku(Sudoku * s)
{
    Sudoku *r = malloc(sizeof(Sudoku));
    r->bdim = s->bdim;
    r->dim = s->dim;
    r->peers_size = s->peers_size;
    r->solved = s->solved;

    r->sol_count = s->sol_count;
    r->unit_list = s->unit_list;
    r->peers = s->peers;

    r->values = calloc(s->dim * s->dim, sizeof(unsigned long long));
    assert(s->values);

    for(int i = 0; i < s->dim * s->dim; i++)
            r->values[i] = s->values[i];

    return r;
}


int fazerTarefas (Sudoku *s, int inicio, int fim, int minI, int possibilidades []) {
    
    Sudoku* vetoresSudoku[fim - inicio];

    #pragma omp parallel
    {

        #pragma omp single nowait
        for (int i = inicio, a = 0; i < fim; i++, a++) {
            vetoresSudoku[a] = copy_sudoku(s);
            vetoresSudoku[a]->status = assign(vetoresSudoku[a], minI, possibilidades[i]);
            if (vetoresSudoku[a]->status) {
                int min = INT_MAX, minI = 0;
                #pragma omp task private(min, minI)
                {
                    encontrarSquareMenosPossibilidades(vetoresSudoku[a], &min, &minI);
                    for (int k = 1; k <= vetoresSudoku[a]->dim; k++) {
                        if (cell_v_get(&vetoresSudoku[a]->values[minI], k)) {
                            search(copy_sudoku(vetoresSudoku[a]), minI, k);
                        }
                    }
                }
            }
        }

        #pragma omp taskwait
    }

    return terminouLocalmente ? 1 : 0;
}

static int search (Sudoku *s, int argMinI, int argK) {
    int i, k;
    if (alguemTerminou || terminouLocalmente) return 0;

    int status = argK > 0 ? assign(s, argMinI, argK) : 1;
    if (!status) return 0;

    int solved = 1;
    for (i = 0; solved && i < s->dim * s->dim; i++) 
        if (cell_v_count(&s->values[i]) != 1) {
            solved = 0;
            break;
        }

    if (solved) {
        s->sol_count++;
        s->status = 1;
        apresentarResultados(s);
        return 1;
    }

    //ok, there is still some work to be done
    int min = INT_MAX;
    int minI = -1;
    int ret = 0;
    
    unsigned long long *values_bkp = malloc (sizeof (unsigned long long) * s->dim * s->dim);

    encontrarSquareMenosPossibilidades(s, &min, &minI);

    if (!jaDividiuProcessos) {
        jaDividiuProcessos = true;

        int qtd_possi = cell_v_count(&s->values[minI]), i, j = 0;
        int possibilidades[qtd_possi];

        for(i = 1 ; i <= s->dim; i++)
        {
            if(cell_v_get(&s->values[minI], i))
            {
                possibilidades[j++] = i;
            }
        }

        int inicio = (int)(world_rank*((double)qtd_possi)/world_size),
            fim = (int)((world_rank + 1)*((double)qtd_possi)/world_size);

        pthread_create(&id, NULL, initPolling, NULL);
        if (fazerTarefas (s, inicio, fim, minI, possibilidades)) return 1;
        pthread_join(id, NULL);
        return 0;

    } else {            
        for (k = 1; k <= s->dim; k++) {
            if (cell_v_get(&s->values[minI], k))  {
                for (i = 0; i < s->dim * s->dim; i++)
                    values_bkp[i] = s->values[i];
                
                if (search (s, minI, k)) {
                    ret = 1;
                    goto FR_RT;
                } else {
                    for (i = 0; i < s->dim * s->dim; i++) 
                        s->values[i] = values_bkp[i];
                }
            }
        }
    }

    FR_RT:
    free(values_bkp);
    
    return ret;
}

int solve(Sudoku *s) {
    return search(s, -1, 0);
}

int main (int argc, char **argv) {

    MPI_Init(NULL, NULL);

    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    MPI_Irecv(&alguemTerminou, 1, MPI_INT, world_rank == 1 ? 0 : 1, 123, MPI_COMM_WORLD, &polling);

    Sudoku *s;

    if (world_rank == 0) {
        int size;
        assert(scanf("%d", &size) == 1);
        assert (size <= MAX_BDIM);
        int buf_size = size * size * size * size;
        unsigned char buf[buf_size];

        for (int i = 0; i < buf_size; i++) {
            if (scanf("%hhu", &buf[i]) != 1) {
                printf("error reading file (%d)\n", i);
                exit(1);
            }
        }

        MPI_Send(&size, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
        MPI_Send(&buf, buf_size, MPI_UNSIGNED_CHAR, 1, 0, MPI_COMM_WORLD);

        s = create_sudoku(size, buf);
    } else {
        int size;

        MPI_Recv(&size, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        int buf_size = size * size * size * size;
        unsigned char buf[buf_size];
        MPI_Recv(&buf, buf_size, MPI_UNSIGNED_CHAR, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);        

        s = create_sudoku(size, buf);
    } 
    if (s) {
        int result = solve(s);
        if (!result && alguemTerminou == 0 && terminouLocalmente == false) printf("Could not solve puzzle.\n");
        destroy_sudoku(s);
    } else {
        printf("Could not load puzzle.\n");
    }

    MPI_Finalize();

    return 0;
}