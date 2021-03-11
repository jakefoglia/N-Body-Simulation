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

#define N 50
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

// global arrays and vars
Particle* PAoS;
Node* TAoS;

uint32_t TAoS_current_size;
uint32_t TAoS_max_size;  

float sys_boundary;
float tree_boundary; 



std::string particleToString(Particle& p)
{
    return std::string("mass : ") + std::to_string(p.mass) + "\npos: (" + std::to_string(p.pos.x) + ", " +  std::to_string(p.pos.y) + ", " +  std::to_string(p.pos.z) + ")"
        + "\nvel: (" +std::to_string(p.vel.x) + ", " +  std::to_string(p.vel.y) + ", " +  std::to_string(p.vel.z) + ")\n\n";
}


int depth(int x) // N = 4 for quadtree, N = 8 for octree
{
    return  floor (log( (7) * x + 1) / log(8) );
}


bool is_child_tree_index(uint32_t pti, uint32_t cti)
{
    return ((cti - 1) / 8 == pti);
}

bool has_children(Node*& parent) // is this equivalent to (!is_leaf) ? 
{
    uint32_t fci = parent->first_child_TAoS_index;

    if(fci > TAoS_current_size)
    {
        printf("SOMETHING IS WRONG!!\n");
        return false; 
    }

    Node* child = &TAoS[fci];
    if(child->tree_index != UINT32_MAX && is_child_tree_index(parent->tree_index, child->tree_index))
        return true;

    return false; 

}
/*
bounds get_bounds(uint32_t tree_index)
{
    
}
*/
bool sub_region_occupied(Node*& parent, float3& pos)
{
    uint32_t ti = parent->tree_index; 
    
    float3 lb;
    float3 ub;
}

uint8_t get_subregion_index(bounds& parent_region, float3& pos) // returns 0 through 7, specifying the subregion
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
    uint8_t index = get_subregion_index(parent_region, pos);
    return get_subregion(parent_region, index);
}



void insert_node(Node*& TAoS, Particle* p)
{
    bounds b = bounds(-tree_boundary, tree_boundary, -tree_boundary, tree_boundary, - tree_boundary, tree_boundary); // start with entire tree region
    Node* current_parent = TAoS;
    while(sub_region_occupied (current_parent, p->pos) )
    {
        
    }
}


void generate_TAoS(Node*& TAoS, Particle*& PAoS)
{
    // initialize root node
    Node* root = (Node*) malloc(sizeof(Node));
    root->tree_index = 0;
    root->PAoS_index = UINT32_MAX; // not a particle
    root->parent_TAoS_index = UINT32_MAX; // no parent
    root->first_child_TAoS_index = 1;
    root->com = float3{0.0f, 0.0f, 0.0f};

    TAoS[0] = *root;
    TAoS_current_size = 1;


    // fill rest of array with 'empty nodes' 
    for(int i = 1; i < TAoS_max_size; i++)
    {
        Node* n = (Node*) malloc(sizeof(Node));
        n->tree_index = n->PAoS_index = UINT32_MAX;
        TAoS[i] = *n;
    }

    for(int i = 0; i < N; i++)
    {
        insert_node(TAoS, &PAoS[i]);
    }
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
        printf(particleToString(*p).c_str() ,"\n");
    }
    sys_com = sys_com * (1.0f/sys_mass);

    sys_boundary = max_dim;
    tree_boundary = max_dim * 1.25f; // 1/2 of the cubic width

    uint32_t fos = 10; // max 9:1 ratio of internal to leaf nodes

    TAoS_max_size = fos * N;
    TAoS = (Node*) malloc (TAoS_max_size * sizeof(Node));

    generate_TAoS(TAoS, PAoS);

}