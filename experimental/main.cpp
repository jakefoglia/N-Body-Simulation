#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <vector> 
#include <time.h>
#include <sstream>
#include <iomanip>
#include <string>
#include <cmath>

//const uint32_t N = 50;
//const float G = 6.6743e-11; 

#define N 1000
#define G 6.6743e-11



typedef struct float3
{
    float3() {};
    float3(float x, float y, float z) : x(x), y(y), z(z) {}

    float x, y, z;
    inline float3 operator*(float s) const { return float3(x*s, y*s, z*s); }
    inline float3 operator+(const float3& a) const { return float3(x+a.x, y+a.y, z+a.z); }

} float3;

typedef struct bounds
{
    bounds() {};
    bounds(float x_lb, float x_ub, float y_lb, float y_ub, float z_lb, float z_ub) : x_lb(x_lb), x_ub(x_ub), y_lb(y_lb), y_ub(y_ub), z_lb(z_lb), z_ub(z_ub) {}
    float x_lb, x_ub, y_lb, y_ub, z_lb, z_ub;

    inline bounds operator*(float s) const { return bounds(x_lb*s, x_ub*s, y_lb*s, y_ub*s, z_lb*s, z_ub*s); }
} bounds;

typedef struct Particle
{
    float3 pos;
    float3 vel;
    float mass; 

    uint32_t TAoS_index; 
} Particle;

typedef struct Node
{
    uint32_t tree_index; 
    uint32_t PAoS_index; 

    uint32_t parent_TAoS_index;
    uint32_t first_child_TAoS_index; 

    float3 com;
    float total_mass; 
} Node;

typedef struct Index_Pair
{
    uint32_t index1;
    uint32_t index2;
} Index_Pair;



// global arrays and vars
Particle* PAoS;
Node* TAoS;
Node* TAoS_swap_buffer;

uint32_t TAoS_current_size;
uint32_t TAoS_max_size;  

float sys_boundary;
float tree_boundary; 


Index_Pair* sort_array;
uint32_t* map;


bool is_internal_node(Node* node);
bool is_particle_node(Node* node);
bool is_internal_node(Node* node);
bool is_blank_node(Node* node);

std::string particleToString(Particle& p)
{
    return std::string("mass : ") + std::to_string(p.mass) + "\npos: (" + std::to_string(p.pos.x) + ", " +  std::to_string(p.pos.y) + ", " +  std::to_string(p.pos.z) + ")"
        + "\nvel: (" +std::to_string(p.vel.x) + ", " +  std::to_string(p.vel.y) + ", " +  std::to_string(p.vel.z) + ")\n\n";
}
std::string nodeToString(Node& n)
{
    return std::string("tree index : ") + std::to_string(n.tree_index) + "\ncom: (" + std::to_string(n.com.x) + ", " +  std::to_string(n.com.y) + ", " +  std::to_string(n.com.z) + ")"
        + "\ntotal mass: " +std::to_string(n.total_mass) + "\n is_particle : " +  std::to_string(is_particle_node(&n)) + "\n\n";
}

int depth(int x) // N = 4 for quadtree, N = 8 for octree
{
    return  floor (log( (7) * x + 1) / log(8) );
}

/* // dont use this
Node* sub_region_occupied(Node* parent, bounds& parent_region, float3& pos)
{
    uint8_t sub_index = get_subregion_index(parent_region, pos);
    uint32_t sub_tree_index = parent->tree_index * 8 + 1 + sub_index; // + 1 bc sub_index is 0 indexed
    uint32_t first_child_index = parent->first_child_TAoS_index;

    for(int i = first_child_index; i < first_child_index + 8 && i < TAoS_current_size; i++)
    {
        if(TAoS[i].tree_index == sub_tree_index)
        {
            return &TAoS[i];
        }
    } 
    return nullptr;
}
*/

/*
bounds get_bounds(uint32_t tree_index)
{
    
}
*/


uint8_t get_subregion_index(bounds& parent_region, float3& pos) // returns 0 through 7, specifying the subregion. x is 3rd bit, y is 2nd bit, z is 1st bit
{
    float3 center;
    center.x = (parent_region.x_lb + parent_region.x_ub) / 2.0f;
    center.y = (parent_region.y_lb + parent_region.y_ub) / 2.0f;
    center.z = (parent_region.z_lb + parent_region.z_ub) / 2.0f;

    uint32_t index = 0;
    if(pos.x >= center.x)
        index = index | (1 << 2);
    if(pos.y >= center.y)
        index = index | (1 << 1);
    if(pos.z >= center.z)
        index = index | (1 << 0);
    
    return index; 
}
bounds get_subregion(bounds& parent_region, uint8_t index)
{
     float3 center;
    center.x = (parent_region.x_lb + parent_region.x_ub) / 2.0f;
    center.y = (parent_region.y_lb + parent_region.y_ub) / 2.0f;
    center.z = (parent_region.z_lb + parent_region.z_ub) / 2.0f;
    bounds region;
    if(index & (1 << 2)) // x upper half
    {
        region.x_lb = center.x;
        region.x_ub = parent_region.x_ub;
    }
    else
    {
        region.x_lb = parent_region.x_lb;
        region.x_ub = center.x;
    }
    if(index & (1 << 1)) // y upper half
    {
        region.y_lb = center.y;
        region.y_ub = parent_region.y_ub;       
    }
    else
    {
        region.y_lb = parent_region.y_lb;
        region.y_ub = center.y;
    }
    if(index & (1 << 2)) // z upper half
    {
        region.z_lb = center.z;
        region.z_ub = parent_region.z_ub;
    }
    else
    {
        region.z_lb = parent_region.z_lb;
        region.z_ub = center.z;
    }
    return region; 
}
bounds get_subregion(bounds& parent_region, float3& pos)
{
    //uint8_t index = get_subregion_index(parent_region, pos);
    //return get_subregion(parent_region, index);

    float3 center;
    center.x = (parent_region.x_lb + parent_region.x_ub) / 2.0f;
    center.y = (parent_region.y_lb + parent_region.y_ub) / 2.0f;
    center.z = (parent_region.z_lb + parent_region.z_ub) / 2.0f;

    bounds region;
    if(pos.x >= center.x)
    {
        region.x_lb = center.x;
        region.x_ub = parent_region.x_ub;
    }
    else
    {       
        region.x_lb = parent_region.x_lb;
        region.x_ub = center.x;
    }
    if(pos.y >= center.y)
    {
        region.y_lb = center.y;
        region.y_ub = parent_region.y_ub;    
    }
    else
    {
        region.y_lb = parent_region.y_lb;
        region.y_ub = center.y;
    }
    if(pos.z >= center.z)
    {
        region.z_lb = center.z;
        region.z_ub = parent_region.z_ub;
    }
    else
    {
        region.z_lb = parent_region.z_lb;
        region.z_ub = center.z;
    }
    return region; 


}

bool is_child_tree_index(uint32_t pti, uint32_t cti)
{
    return ((cti - 1) / 8 == pti);
}

/* method is broken and useless. Do not use
bool has_children(Node*& parent) // do blanks count as chilren? 
{
    uint32_t fci = parent->first_child_TAoS_index;

    if(fci == UINT32_MAX)
        return false; 

    if(fci > TAoS_current_size)
    {
        printf("SOMETHING IS WRONG!!\n");
        return false; 
    }

    Node* child = &TAoS[fci];
    


    // keep this for safety for now
    if(child->tree_index != UINT32_MAX && !is_child_tree_index(parent->tree_index, child->tree_index))
    {
        printf("\nSOMETHING IS WRONG WITH TAoS child indexing!!!\n");
        printf("fci %i\n", fci);
        printf("pti %i\n", parent->tree_index);
        printf("cti %i\n", child->tree_index);
        
        printf("is child : %i\n", ((child->tree_index - 1) / 8 == parent->tree_index));
    }

    if(child->tree_index != UINT32_MAX && is_child_tree_index(parent->tree_index, child->tree_index))
        return true;

    return false; 

}
*/

bool is_internal_node(Node* node)
{
            // cant be empty node             // cant be particle               // must have children
    return(node->tree_index != UINT32_MAX && node->PAoS_index == UINT32_MAX && node->first_child_TAoS_index != UINT32_MAX);
}
bool is_blank_node(Node* node) // a blank is a leaf node thats not a particle. 
{
            // cant be empty node             // cant be particle               // cant have children
    return (node->tree_index != UINT32_MAX && node->PAoS_index == UINT32_MAX && node->first_child_TAoS_index == UINT32_MAX);
}

bool is_particle_node(Node* node)
{       
            // cant be empty node           // must have index on the PAoS array   // cant have children
    return (node->tree_index != UINT32_MAX && node->PAoS_index != UINT32_MAX    && node->first_child_TAoS_index == UINT32_MAX);
}
bool is_empty_node(Node* node)
{
            // all nonempty nodes must have tree index
    return node->tree_index == UINT32_MAX; 
}


// 3 types of nodes: internal, particle, blank. 
// root starts as blank: it has no children but its not a particle
// then after trying to add one particle, we add the particle and 7 other blank node children which will be turned into particles then internal nodes when necessary. 



void insert_blanks(uint32_t parent_TAoS_index)
{
    Node* parent = &TAoS[parent_TAoS_index];
    parent->first_child_TAoS_index = TAoS_current_size;

    // insert 8 blank children
    for(int i = 1; i <= 8; i++)
    {
        Node* n = &TAoS[TAoS_current_size++];

        n->PAoS_index = UINT32_MAX;
        n->tree_index = parent->tree_index*8 + i;
        n->first_child_TAoS_index = UINT32_MAX; // no children yet, these are blanks
        n->parent_TAoS_index = parent_TAoS_index;

        n->com = float3{0, 0, 0};
        n->total_mass = {0};
    }
    //printf("inserted blanks as child of TAoS[%i]\n", parent_TAoS_index); // debug
}


void insert_particle_on_blank(Node* blank,  uint32_t PAoS_index) //uint32_t parent_TAoS_index,
{
    if(!is_blank_node(blank))
    {
        printf("error in insert_particle_on_blank : node isn't blank!!!\n");
        return;
    }

    Particle* p =  &PAoS[PAoS_index];

    blank->com = p->pos;
    blank->total_mass = p->mass;
    blank->PAoS_index = PAoS_index; 
    blank->first_child_TAoS_index = UINT32_MAX;
    //blank->parent_TAoS_index = parent_TAoS_index; // is this necessary? when blanks were inserted this should have been done already
}

void insert_internal_on_particle(uint32_t TAoS_index, bounds& b)
{
    Node bkp = TAoS[TAoS_index]; // create a copy
    //uint32_t parent_TAoS_index = bkp.parent_TAoS_index;

    insert_blanks(TAoS_index);

    uint8_t sub_index = get_subregion_index(b, bkp.com);
    uint32_t target_TAoS_index = TAoS[TAoS_index].first_child_TAoS_index + sub_index; 

    insert_particle_on_blank(&TAoS[target_TAoS_index], bkp.PAoS_index); // theres gotta be a better way than this. override the method to just take a Node& particle 

}




void insert_particle_in_tree(Node* TAoS, uint32_t PAoS_index) // a node should have 0 children, or 8 children (must be consecutive)
{
    Particle* p = &PAoS[PAoS_index];
    bounds current_bounds = bounds(-tree_boundary, tree_boundary, -tree_boundary, tree_boundary, -tree_boundary, tree_boundary); // start with entire tree boundary
    Node* current_node = TAoS;
    Node* parent_node = TAoS; 

    //uint32_t parent_TAoS_index;

    bool flag = false; 
    
    /*
    // traverse until we find a blank or a particle to split
    while(!is_blank_node(current_node) && !is_particle_node(current_node) ) // then it should have 8 consecutive children
    */

    // loop until we find a blank. branch tree if we hit a particle. 
    while(!is_blank_node(current_node))  // while(is_internal) ? is that equivalent
    {   
        // TODO: add COM and total_mass contribution on the way down the tree
        //current_node->com = ( (current_node->com * current_node->total_mass) + (p->mass * p->pos) ) / (current_node->total_mass + p->mass); 

        // advance the tree
        if(flag)
        {
            parent_node = current_node;    
        }
        flag = true;

        uint8_t sub_index = get_subregion_index(current_bounds, p->pos);
        current_bounds =  get_subregion(current_bounds, sub_index);

        
        if( (!is_particle_node(current_node) && current_node->first_child_TAoS_index == UINT32_MAX) || current_node->first_child_TAoS_index == 0) // debug
            printf("uhh ohh in insert_node()\n\n");


        if(is_particle_node(current_node))
        {
            insert_internal_on_particle(current_node - TAoS, current_bounds); // ?????? probably not right
        }
        else // only advance if it was already and internal node. 
        {
            // this is why we must insert 8 children at once, or none 
            current_node = &TAoS[current_node->first_child_TAoS_index + sub_index]; // is this causing overflow ? 
        }

    }

    uint32_t parent_TAoS_index = (parent_node - TAoS); // is this right?
    if(is_blank_node(current_node))
    {
        
        //printf("inserting PAoS[%i] on blank child of tree_index %i\n", PAoS_index, TAoS[parent_TAoS_index].tree_index); // debug
        insert_particle_on_blank(current_node, PAoS_index);
        return; 
    }

 



/*
    Node* occupying_node = sub_region_occupied (current_node, current_bounds, p->pos);
    while(occupying_node != nullptr)
    {
        bool isParticle = (occupying_node->PAoS_index != UINT32_MAX); 
        if(isParticle)
        {
            // branch the tree
            break;
        }
        else
        {
            // TODO: add COM and total_mass contribution on the way down the tree
            //current_node->com = ( (current_node->com * current_node->total_mass) + (p->mass * p->pos) ) / (current_node->total_mass + p->mass); 

            // advance
            current_bounds =  get_subregion(current_bounds, p->pos);
            current_node = occupying_node;

            
            
        }
    } 

    */
}


void generate_TAoS(Node* TAoS, Particle* PAoS)
{
    TAoS_current_size = 0;

    // initialize root node
    Node* root = TAoS;
    root->tree_index = 0;
    root->PAoS_index = UINT32_MAX; // not a particle
    root->parent_TAoS_index = UINT32_MAX; // no parent
    root->first_child_TAoS_index = 1; // we will put 8 blanks there
    root->com = float3{0.0f, 0.0f, 0.0f};

    
    TAoS_current_size++;

    insert_blanks(0);

/*
    // insert 8 blank children
    for(int ti = 1; ti <= 8 ; ti++)
    {
        Node* n = (Node*) malloc(sizeof(Node));
        n->PAoS_index = UINT32_MAX;
        n->tree_index = ti;
        n->first_child_TAoS_index = UINT32_MAX; // no children yet, these are blanks
        n->parent_TAoS_index = 0;
        TAoS[ti] = *n;

        TAoS_current_size++;
    }
*/
    // fill rest of array with 'empty nodes' 
    for(int ti = 9; ti < TAoS_max_size; ti++)
    {
        Node* n = &TAoS[ti];
        n->tree_index = UINT32_MAX;
        n->PAoS_index = UINT32_MAX;
        n->first_child_TAoS_index = UINT32_MAX;
        n->parent_TAoS_index = UINT32_MAX;
    }

    for(int pi = 0; pi < N; pi++)
    {
        insert_particle_in_tree(TAoS, pi);
        //printf("inserted p\n\n");
    }
}

// TODO: comapre with std::sort on the TAoS. how will we get the map to update the parent and child indices? 

void sort_TAoS(Node*& TAoS, Node*& TAoS_swap_buffer, Index_Pair* n2o, uint32_t* m, uint32_t n, uint32_t max) 
{
    for(uint32_t i = 0; i < n; i++)
    {
        n2o[i].index1=TAoS[i].tree_index;
        n2o[i].index2=i; // TAoS_index_old
    }
    
    // insertion sort
    Index_Pair key; 
    uint32_t i, j; 
    for (i = 1; i < n; i++)
    { 
        key = n2o[i]; 
        j = i - 1; 
 
        while (j >= 0 && n2o[j].index1 > key.index1)
        { 
            n2o[j + 1] = n2o[j]; 
            j = j - 1; 
        } 
        n2o[j + 1] = key; 
    }
    
    // fill map based on inverse  (can we do this inside the insertion sort ? )
    for(uint32_t i = 0; i < n; i++)
    {
        m[n2o[i].index2] = i; // map [TAoS_index_old] = TAoS_index_new
    }

    // root is special case ( has no map for parent_TAoS_index)
    TAoS_swap_buffer[0] = TAoS[0];
    TAoS_swap_buffer[0].first_child_TAoS_index = m[TAoS[0].first_child_TAoS_index];

    // remap remaining non-empty nodes
    for(int i = 1; i < n; i++)
    {
        TAoS_swap_buffer[m[i]] = TAoS[i];
        if(TAoS[i].first_child_TAoS_index < UINT32_MAX)
            TAoS_swap_buffer[m[i]].first_child_TAoS_index = m[TAoS[i].first_child_TAoS_index]; // if this seg faults, we are trying to map based on an invalid index

        if(TAoS[i].parent_TAoS_index < UINT32_MAX)
            TAoS_swap_buffer[m[i]].parent_TAoS_index = m[TAoS[i].parent_TAoS_index]; // if this seg faults, we are trying to map based on an invalid index
    }
    for(int i = n; i < TAoS_max_size; i++)
    {
        TAoS_swap_buffer[i] = TAoS[i];
    }

    // swap the buffers
    Node* temp = TAoS;
    TAoS = TAoS_swap_buffer;
    TAoS_swap_buffer = temp;
}

int main ()
{
    float3 sys_com = float3{0.0f, 0.0f, 0.0f};
    float sys_mass = 0.0f;
    float max_dim = 0.0f;

    //std::vector<Particle*> particles;
    PAoS = (Particle*) malloc(N * sizeof(Particle)); 

    srand(time(NULL));
    uint32_t range = 1000; 
    for(int i = 0; i < N; i++)
    {
        Particle* p = (Particle*) malloc(sizeof(Particle)) ;
        p->pos=float3( rand() % (range) - range/2.0f, rand() % range - range/2.0f, rand() % range - range/2.0f);
        p->vel=float3(0,0,0);
        p->mass = 1;
        
        if(std::abs(p->pos.x) > max_dim)
            max_dim = p->pos.x;
        
        if(std::abs(p->pos.y) > max_dim)
            max_dim = p->pos.y;
        
        if(std::abs(p->pos.z) > max_dim)
            max_dim = p->pos.z;

        sys_com = sys_com + ((p->pos) * p->mass );
        sys_mass += p->mass;

        //printf("p_pos %f %f %f\n", p->pos.x, p->pos.y, p->pos.z);
        //printf("s_com %f %f %f\n", sys_com.x, sys_com.y, sys_com.z);

        PAoS[i] = *p;
        //printf(particleToString(*p).c_str() ,"\n");
    }
    sys_com = sys_com * (1.0f/sys_mass);

    sys_boundary = max_dim;
    tree_boundary = max_dim * 1.25f; // 1/2 of the cubic width

    uint32_t fos = 10; // max 9:1 ratio of internal to leaf nodes

    TAoS_max_size = fos * N;
    TAoS = (Node*) malloc (TAoS_max_size * sizeof(Node));
    TAoS_swap_buffer = (Node*) malloc (TAoS_max_size * sizeof(Node));
    

    Node* t0 = TAoS;
    Node* t1 = &TAoS[1];
    Node* t2 = &TAoS[2];
    Node* t3 = &TAoS[3];
    Node* t4 = &TAoS[4];
    Node* t5 = &TAoS[5];
    Node* t6 = &TAoS[6];
    Node* t7 = &TAoS[7];
    Node* t8 = &TAoS[8];

    sort_array = (Index_Pair*) malloc(TAoS_max_size * sizeof(Index_Pair));
    map = (uint32_t*) malloc(TAoS_max_size * sizeof(uint32_t));


    generate_TAoS(TAoS, PAoS);

    sort_TAoS(TAoS, TAoS_swap_buffer, sort_array, map, TAoS_current_size, TAoS_max_size);

    for(int i = 0; i < TAoS_current_size + 1; i++)
    {
        printf("ti: %i\n", TAoS[i].tree_index);
        //printf(nodeToString(TAoS[i]).c_str());
    }


    printf("\n");
    printf("PAoS size :  %i\n", N);
    printf("TAoS size :  %i\n", TAoS_current_size);
    printf("ratio     :  %4.2f\n", ((TAoS_current_size * 1.0f) / N));


    /*
   bounds b1{-1.0f, 1.0f, -1.0f, 1.0f, -1.0f, 1.0f};
   float3 pos = {-0.5f, 0.5f, 0.5f};

   printf("sub_index : %i\n\n", get_subregion_index(b1, pos));


   printf("index: %i ", (&TAoS[26] - TAoS) );
*/
}

/*
Last Error:

Jake@JakeASUS MINGW64 /d/Git/N-Body-Simulation/experimental (master)
$ ./a.exe
PAoS size :  50
TAoS size :  273
ratio     :  5.46

Jake@JakeASUS MINGW64 /d/Git/N-Body-Simulation/experimental (master)
$ ./a.exe
error in insert_particle_on_blank : node isn't blank!!!
uhh ohh in insert_node()

Segmentation fault


*/



/*

PLAN FOR SORTING TAoS

Generate a new AoS. Each struct holds this info from a node on the TAoS:

Sorting Struct: 
{
    tree_index  (for sorting)
    TAoS_index 
}

We radix sort based on tree index. 

This gives us a map 
    from    TAoS_index_new    (location in sorted array)
    to      TAoS_index_old    (TAoS_index entry in the struct)

We want the inverse of this map.
So we create a new array for the inverse map

the location in the array will represent the TAoS_index_old
the values will hold the TAoS_index_new

So for each entry in the sorted array, we read the TAoS_index (the old index locations)
and we use this index to go to that location in the new map.

Then at that location, we write the index in the sorted array of where we found it.

This will give us our desried map
    from    TAoS_index_old    (location in map)
    to      TAoS_index_new    (value at that location)


Using this map we can easily sort our TAoS efficiently.
Then after sorting, we go to each node and update their:
    TAoS_parent_index
 &  TAoS_first_child_index 
 
 based on the map   

*/

/*
PDF paper format:
YYYY-where published-title
*/
