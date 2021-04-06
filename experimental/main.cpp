#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <vector> 
#include <sstream>
#include <iomanip>
#include <string>
#include <cmath>
#include <stack>
#include <ctime>
#include <chrono>
#include <pthread.h> 
#include <assert.h>

//const uint32_t N = 50;
//const float G = 6.6743e-11; 

#define N 1000
#define G 1.0f
#define dt  0.1f
#define theta  1.0f
#define epsilon 0.01f
#define fos 20
#define boundary_fos 1.25f

typedef struct float3
{
    float3() {};
    float3(float x, float y, float z) : x(x), y(y), z(z) {}

    float x, y, z;
    inline float3 operator*(float s) const { return float3(x*s, y*s, z*s); }
    inline float3 operator/(float s) const { return float3(x/s, y/s, z/s); }
    inline float3 operator+(const float3& a) const { return float3(x+a.x, y+a.y, z+a.z); }

} float3;

typedef struct bounds
{
    bounds() {};
    bounds(float x_lb, float x_ub, float y_lb, float y_ub, float z_lb, float z_ub) : x_lb(x_lb), x_ub(x_ub), y_lb(y_lb), y_ub(y_ub), z_lb(z_lb), z_ub(z_ub) {}
    float x_lb, x_ub, y_lb, y_ub, z_lb, z_ub;

    inline bounds operator*(float s) const { return bounds(x_lb*s, x_ub*s, y_lb*s, y_ub*s, z_lb*s, z_ub*s);  }
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

typedef struct xthread_arg_struct 
{
    Node* TAoS;
    Particle* PAoS_start_address;
    uint32_t n;
} xthread_arg_struct;



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
bool is_blank_node(Node* node);
bool is_empty_node(Node* node);


std::string particleToString(Particle& p)
{
    return std::string("mass : ") + std::to_string(p.mass) + "\npos: (" + std::to_string(p.pos.x) + ", " +  std::to_string(p.pos.y) + ", " +  std::to_string(p.pos.z) + ")"
        + "\nvel: (" +std::to_string(p.vel.x) + ", " +  std::to_string(p.vel.y) + ", " +  std::to_string(p.vel.z) + ")\n\n";
}
std::string nodeToString(Node& n)
{
    std::string type ="";
    if(is_empty_node(&n))
        type+="empty ";
    if(is_blank_node(&n))
        type+="blank ";
    if(is_internal_node(&n))
        type+="internal ";
    if(is_particle_node(&n))
        type+="particle ";    

    return std::string("tree index : ") + std::to_string(n.tree_index) + "  " + type + "\ncom: (" + std::to_string(n.com.x) + ", " +  std::to_string(n.com.y) + ", " +  std::to_string(n.com.z) + ")"
        + "\ntotal mass: " +std::to_string(n.total_mass) +  "\n\n";
}

int depth(int x) // N = 4 for quadtree, N = 8 for octree
{
    return  floor (log( (7) * x + 1) / log(8) );
}


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
/*
bool is_valid(Node* node)
{
    return (node - TAoS < TAoS_current_size);
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
    return (node->tree_index == UINT32_MAX); // index out of bounds means we treat as empty
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

        n->total_mass = 0.0f;
        n->com = float3{0.0f, 0.0f, 0.0f};
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
    Node* node = &TAoS[TAoS_index];
    if(!is_particle_node(node))
    {
        printf("error in insert_inernal_on_particle : node isn't particle!!!\n");
        return;
    }

    Node bkp = *node; // create a copy
    //uint32_t parent_TAoS_index = bkp.parent_TAoS_index;

    insert_blanks(TAoS_index);
    node->com = bkp.com;
    node->total_mass = bkp.total_mass;
    node->PAoS_index = UINT32_MAX;

    uint8_t sub_index = get_subregion_index(b, bkp.com);
    uint32_t target_TAoS_index = node->first_child_TAoS_index + sub_index; 

    insert_particle_on_blank(&TAoS[target_TAoS_index], bkp.PAoS_index); // theres gotta be a better way than this. override the method to just take a Node& particle 

}

void merge_particles(Node* existing_particle,  uint32_t PAoS_index)
{
     if(!is_particle_node(existing_particle))
    {
        printf("error in merge_particles : existing node isn't a particle!!!\n");
        return;
    }

    Particle* p =  &PAoS[PAoS_index];

    existing_particle->com = ((existing_particle->com * existing_particle->total_mass) + p->pos * p->mass) / (existing_particle->total_mass + p->mass);
    existing_particle->total_mass = (existing_particle->total_mass + p->mass);
    printf("merged new mass %6.4f\n", existing_particle->total_mass);

    //existing_particle->PAoS_index = PAoS_index; // can we just leave the PAoS index the same as the original  ? 

}



void insert_particle_in_tree(Node* TAoS, uint32_t PAoS_index) // a node should have 0 children, or 8 children (must be consecutive)
{
 
    //printf("\ninserting PAoS[%i] TCS: %i  TMS: %i  POS(%6.4f, %6.4f, %6.4f)\n", PAoS_index, TAoS_current_size, TAoS_max_size, PAoS[PAoS_index].pos.x, PAoS[PAoS_index].pos.y, PAoS[PAoS_index].pos.z);
    
    Particle* p = &PAoS[PAoS_index];
    bounds current_bounds = bounds(-tree_boundary, tree_boundary, -tree_boundary, tree_boundary, -tree_boundary, tree_boundary); // start with entire tree boundary
    Node* current_node = TAoS;
    Node* parent_node = TAoS; 

    bool oob = false;
    if(PAoS[PAoS_index].pos.x < current_bounds.x_lb || PAoS[PAoS_index].pos.x > current_bounds.x_ub ||
    PAoS[PAoS_index].pos.y < current_bounds.y_lb || PAoS[PAoS_index].pos.y > current_bounds.y_ub ||
    PAoS[PAoS_index].pos.z < current_bounds.z_lb || PAoS[PAoS_index].pos.z > current_bounds.z_ub) 
    {
        printf("PARTICLE IS OUT OF BOUNDS!!!!!!!!!\n"); 
        printf("(%6.4f, %6.4f, %6.4f)\t(%6.4f, %6.4f, %6.4f) - (%6.4f, %6.4f, %6.4f)\n", PAoS[PAoS_index].pos.x, PAoS[PAoS_index].pos.y, PAoS[PAoS_index].pos.z, current_bounds.x_lb, current_bounds.y_lb, current_bounds.z_lb, 
                        current_bounds.x_ub, current_bounds.y_ub, current_bounds.z_ub);
        oob = true;
    }


    //uint32_t parent_TAoS_index;

    bool flag = false; 
    
    /*
    // traverse until we find a blank or a particle to split
    while(!is_blank_node(current_node) && !is_particle_node(current_node) ) // then it should have 8 consecutive children
    */

    // loop until we find a blank. branch tree if we hit a particle. 
    //while(!is_blank_node(current_node))  // while(is_internal) ? is that equivalent
    uint32_t depth = 0;
    uint32_t pinarow = 0;
    bool merged = false; 
    while(is_internal_node(current_node) || is_particle_node(current_node) && !merged)
    {   
        // TODO: add COM and total_mass contribution on the way down the tree
        //current_node->com = ( (current_node->com * current_node->total_mass) + (p->mass * p->pos) ) / (current_node->total_mass + p->mass); 

        //printf("TAoS[%i] bounds(%6.4f - %6.4f, %6.4f - %6.4f, %6.4f - %6.4f)\t", current_node - TAoS, current_bounds.x_lb, current_bounds.x_ub, current_bounds.y_lb, current_bounds.y_ub, current_bounds.z_lb, current_bounds.z_ub);
        //printf("TAoS[%i]\t", current_node - TAoS);
        // advance the tree
        //printf(nodeToString(*current_node).c_str());
        
        if(depth != 0)
        {
            parent_node = current_node;
        }
        depth++;
        /*
        if(flag)
        {
            parent_node = current_node;    
        }
        flag = true;
        */
        
        uint8_t sub_index = get_subregion_index(current_bounds, p->pos);
        current_bounds =  get_subregion(current_bounds, sub_index);

        if(current_node->first_child_TAoS_index == 0)
            printf("uhh ohh 1 in insert_node()\n\n");
        if( (!is_particle_node(current_node) && current_node->first_child_TAoS_index == UINT32_MAX) ) // debug
            printf("uhh ohh 2 in insert_node()\n\n");


        if(is_particle_node(current_node))
        {
            pinarow++;
            if(pinarow == 30) // problem is that two particle have very similar positions. just combine them into a single particle node on the tree. 
            {
                
                printf("\npinarow getting big : %i\n", pinarow);

                merge_particles(current_node, PAoS_index);
                merged = true;

                //Node* tmp = TAoS;
                //bounds b = bounds(-tree_boundary, tree_boundary, -tree_boundary, tree_boundary, -tree_boundary, tree_boundary);

                /*
                printf("PAoS[%i] : (%6.4f, %6.4f, %6.4f)\n", PAoS_index, p->pos.x, p->pos.y, p->pos.z);
                Node* tmp = current_node;
                while(tmp != TAoS)
                {
                    printf("com : %6.4f, %6.4f, %6.4f\n", tmp->com.x, tmp->com.y, tmp->com.z);
                    printf("mass: %6.4f\n", tmp->total_mass);
                    tmp = &TAoS[tmp->parent_TAoS_index];
                }
                */
            }
            else
                insert_internal_on_particle(current_node - TAoS, current_bounds); 
        }
        else // only advance if it was already and internal node. 
        {
            // add mass contribution
            current_node->com = ((current_node->com * current_node->total_mass) + (p->pos * p->mass) ) *  (1.0f / (current_node->total_mass + p->mass));
            current_node->total_mass += p->mass; 

            // this is why we must insert 8 children at once, or none 
            current_node = &TAoS[current_node->first_child_TAoS_index + sub_index]; 
        }
    }
    if(oob)
    {
        printf("PAoS[%i] was out of bounds and reached depth %i\n", PAoS_index, depth);
    }

    uint32_t parent_TAoS_index = (parent_node - TAoS); 
    if(is_blank_node(current_node))
    {
        //printf("reached valid blank! Trying to insert at TAoS[%i]\n", current_node-TAoS);
        //printf("inserting PAoS[%i] on blank child of tree_index %i\n", PAoS_index, TAoS[parent_TAoS_index].tree_index); // debug
        insert_particle_on_blank(current_node, PAoS_index);
        
        return; 
    }
    else if (merged)
    {
        printf ("successfully merged particles:\n PAoS[%i] : (%6.4f, %6.4f, %6.4f)\n PAoS[%i] : (%6.4f, %6.4f, %6.4f)\n", current_node->PAoS_index, 
            PAoS[current_node->PAoS_index].pos.x, PAoS[current_node->PAoS_index].pos.y, PAoS[current_node->PAoS_index].pos.z, PAoS_index, PAoS[PAoS_index].pos.x, PAoS[PAoS_index].pos.y, PAoS[PAoS_index].pos.z);
        printf("new total mass: %6.4f\n", current_node->total_mass);
    }
    
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
    root->total_mass = 0.0f;

    
    TAoS_current_size++;

    // insert 8 blank children
    insert_blanks(0);

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
        //printf("%i %i\n",TAoS_current_size, TAoS_max_size);
        insert_particle_in_tree(TAoS, pi);
    }
    //printf("TAoS_current_size %i\n", TAoS_current_size);
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
float norm(float3 vec)
{
    return std::sqrt(vec.x*vec.x + vec.y*vec.y + vec.z*vec.z);
}
void dfs_traverse(Node* TAoS, Particle* p)
{
    ///printf("\n--------------------------------------------------------------\n");
    ///printf("%i  (%6.4f, %6.4f, %6.4f)\n",p - PAoS, p->pos.x, p->pos.y, p->pos.z);

    float3 accel(0.0f, 0.0f, 0.0f);


    // init stack with root
    std::stack<std::pair<Node*, float>> stack;
    stack.push(std::pair<Node*, float>(TAoS, 2.0f*tree_boundary)); 

    // traverse
    while(!stack.empty())
    {
        std::pair<Node*, float> pair = stack.top();
        stack.pop();

        Node* current = pair.first;
        float width = pair.second;


        float3 r = current->com + p->pos * -1.0f; 
        float r_norm = norm(r);

        ///printf("\t%i\tinternal: %i\t r_norm/width %6.2f\t total_mass %6.2f\t", current->tree_index, is_internal_node(current), r_norm / width, current->total_mass );

        if(is_internal_node(current) && (r_norm  / width) < 1.0f ) 
        {
            ///printf("push children");
            // push valid children to stack
            Node* child;
            for(int i = 0; i < 8; i++)
            {
                child = &TAoS[current->first_child_TAoS_index + i];
                if(!is_empty_node(child) && !is_blank_node(child))
                    stack.push(std::pair<Node*, float> (child, width/2.0f));
            }
        }
        else 
        {   
            ///printf("FORCE "); 
            ///printf((is_particle_node(current)) ? " (particle) " : " (internal) " );
            
    
            float3 a(0.0f, 0.0f, 0.0f);
            a = ((r / (r_norm + epsilon)) * (float) (G * current->total_mass)) / ((r_norm + epsilon ) * (r_norm + epsilon));
            accel =  a + accel;
            
            ///printf("accel += (%10.8f, %10.8f, %10.8f) = (%10.8f, %10.8f, %10.8f)", a.x, a.y, a.z, accel.x, accel.y, accel.z );
        }
        ///printf("\n");
    }

    ///printf("\npos: (%10.8f, %10.8f, %10.8f)\n", p->pos.x, p->pos.y, p->pos.z);
    ///printf("vel: (%10.8f, %10.8f, %10.8f)\n", p->vel.x, p->vel.y, p->vel.z);

    ///printf("\naccel: (%10.8f, %10.8f, %10.8f)\n", accel.x, accel.y, accel.z);

    p->vel = p->vel + accel * dt;
    p->pos = p->pos + p->vel*dt; 
    
    // fixed here
    if(std::abs(p->pos.x) > sys_boundary)
        sys_boundary = std::abs(p->pos.x);
        
    if(std::abs(p->pos.y) > sys_boundary)
        sys_boundary = std::abs(p->pos.y);
        
    if(std::abs(p->pos.z) > sys_boundary)
        sys_boundary = std::abs(p->pos.z);
    
    tree_boundary = sys_boundary * boundary_fos; // only safe if we regenerate tree from scratch


    ///printf("\npos: (%10.8f, %10.8f, %10.8f)\n", p->pos.x, p->pos.y, p->pos.z);
    ///printf("vel: (%10.8f, %10.8f, %10.8f)\n", p->vel.x, p->vel.y, p->vel.z);    
    
}
void calculate_forces(Node* TAoS, Particle* PAoS, uint32_t n)
{
    for(int i = 0; i < n; i++)
    {
        dfs_traverse(TAoS, &PAoS[i]); // traverse starting at the root node
    }
}
void* xthreaded_traversal(void* arg)
{
    xthread_arg_struct* arg_struct = (xthread_arg_struct*) arg;

    Node* TAoS = arg_struct->TAoS;
    Particle* PAoS_start_address = arg_struct->PAoS_start_address;
    uint32_t n = arg_struct->n;

    calculate_forces(TAoS, PAoS_start_address, n);

    pthread_exit(0);
    return nullptr;

}
void xthreaded_calculate_forces(uint8_t num_threads, Node* TAoS, Particle* PAoS, uint32_t n)
{
    xthread_arg_struct arg_structs[num_threads];
    if(num_threads < 2)
    {
        calculate_forces(TAoS, PAoS, n);
    }
    else
    {
        pthread_t tids[num_threads-1];
        Particle* running_PAoS_ptr = PAoS;
        uint32_t ppt = n / (num_threads - 1); // particle per thread
        for(int i = 0; i < num_threads - 1; i++)
        {
            arg_structs[i].TAoS = TAoS;
            arg_structs[i].PAoS_start_address = running_PAoS_ptr;
            arg_structs[i].n = ppt;

            pthread_attr_t attr;
            pthread_attr_init(&attr);
          
            pthread_create(&tids[i], &attr, xthreaded_traversal, (void*) &arg_structs[i]);
            running_PAoS_ptr = &running_PAoS_ptr[ppt];
        }
        calculate_forces(TAoS, running_PAoS_ptr, n - ppt*(num_threads-1)); // traverse the remainder in the meantime;

        for(int i = 0; i < num_threads-1; i++)
        {
            pthread_join(tids[i], NULL);
        }
    }
    return;
}


int main (int argc, char** argv)
{
    uint8_t num_threads;
    if(argc < 2) {
        printf("Usage: %s <num_threads>\n", argv[0]);
        printf("Applying Default: <num_threads> = 1\n\n");
    }
    else {
        num_threads = atoi(argv[1]);
        printf("<num_threads> = %i\n", (int) num_threads);
    }

    uint32_t *TMS, *TCS;
    TMS = &TAoS_max_size;
    TCS = &TAoS_current_size;



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
        //p->pos=float3( rand() % (range) - range/2.0f, rand() % range - range/2.0f, rand() % range - range/2.0f);
        p->pos=float3(-500.0f+10.0f*i + rand()%10, -500.0f+10.0f*i +rand()%10, -500.0f+10.0f*i + rand()%10);
        p->vel=float3(0,0,0);
        p->mass = 1;

        if(i == 0) // manually insert large body
        {
            p->pos=float3(10.0f, 10.0f, 10.0f);
            p->vel = float3(0.0f, 0.0f, 0.0f);
            p->mass = 10000.0f;
        }
        
        if(std::abs(p->pos.x) > max_dim)
            max_dim = std::abs(p->pos.x);
        
        if(std::abs(p->pos.y) > max_dim)
            max_dim = std::abs(p->pos.y);
        
        if(std::abs(p->pos.z) > max_dim)
            max_dim = std::abs(p->pos.z);

        sys_com = sys_com + ((p->pos) * p->mass );
        sys_mass += p->mass;

        //printf("p_pos %f %f %f\n", p->pos.x, p->pos.y, p->pos.z);
        //printf("s_com %f %f %f\n", sys_com.x, sys_com.y, sys_com.z);

        PAoS[i] = *p;
        //printf(particleToString(*p).c_str() ,"\n");
    }
    sys_com = sys_com * (1.0f/sys_mass);

    sys_boundary = max_dim;
    //tree_boundary = max_dim * 1.25f; // 1/2 of the cubic width
    tree_boundary = max_dim * boundary_fos;


    TAoS_max_size = fos * N;
    TAoS = (Node*) malloc (TAoS_max_size * sizeof(Node));
    //TAoS_swap_buffer = (Node*) malloc (TAoS_max_size * sizeof(Node)); // we dont need this unless we are sorting

    sort_array = (Index_Pair*) malloc(TAoS_max_size * sizeof(Index_Pair));
    map = (uint32_t*) malloc(TAoS_max_size * sizeof(uint32_t));



    uint64_t update_TAoS_time = 0;
    uint64_t sort_TAoS_time = 0;
    uint64_t calculate_forces_time = 0;
    uint64_t gen_TAoS_time = 0;
    uint64_t total = 0;


    auto t0 =  std::chrono::steady_clock::now();
    generate_TAoS(TAoS, PAoS);
    gen_TAoS_time = (std::chrono::steady_clock::now() - t0).count();
    
    uint32_t max_TAoS_size_reached = 0;
    int steps = 1000;
    for(int i = 0; i < steps; i++)
    {
        auto t1 = std::chrono::steady_clock::now();
        generate_TAoS(TAoS, PAoS);
        //update_TAoS(TAoS, PAoS);
        auto t2 = std::chrono::steady_clock::now();
        //sort_TAoS(TAoS, TAoS_swap_buffer, sort_array, map, TAoS_current_size, TAoS_max_size); // this is optional for cpu version. It just slows us down unecessarily! 
        auto t3 = std::chrono::steady_clock::now();
        calculate_forces(TAoS, PAoS, N);
        //xthreaded_calculate_forces(num_threads, TAoS, PAoS, N);
        auto t4 = std::chrono::steady_clock::now();


        if(TAoS_current_size > max_TAoS_size_reached)	
            max_TAoS_size_reached = TAoS_current_size;	

        update_TAoS_time += static_cast<uint64_t>((t2 - t1).count());
        sort_TAoS_time += static_cast<uint64_t>((t3 - t2).count());
        calculate_forces_time += static_cast<uint64_t>((t4 - t3).count());
    }
    
    update_TAoS_time /= steps;
    sort_TAoS_time /= steps;
    calculate_forces_time /= steps;
    
    total = update_TAoS_time + sort_TAoS_time + calculate_forces_time;


    printf("\n\n\n");
    printf("PAoS size :  %i\n", N);
    printf("TAoS size :  %i\n", TAoS_current_size);
    printf("ratio     :  %4.2f\n", ((TAoS_current_size * 1.0f) / N));

    printf("\nAvg Times:\n");
    printf("generate_TAoS_time : %lf ms\n", gen_TAoS_time/1E6);
    printf("update_TAoS_time : %lf ms\n", update_TAoS_time/1E6);
    printf("sort_TAoS_time : %lf ms\n", sort_TAoS_time/1E6);
    printf("calculate_forces_time : %lf ms\n", calculate_forces_time/1E6);
    printf("total iteration time : %lf ms\n", total/1E6);

    delete[] TAoS;
    delete[] TAoS_swap_buffer;
    delete[] PAoS;
    delete[] sort_array;
    delete[] map;

    return 1;

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

//TODO: deal with out of bounds



// -O3, microway, check for correctness with pthread, test larger dataset, new design ideas

/*
Idea for handling out of bounds without regenerating tree (dynamic adjustment)

when particle is oob, rather than branching infinitely when another oob particle competes with it, we just maintain a linked list for those particles.
we cam't really do the hashing strategy bc then we will unjustly occupy  positions for valid particles (in bounds).

oob linked list! or array ? 
maybe we double tree boundary if the array gets filled, in addition to sys boundary doubling or halving 


remap of tree indices is easier than you think. 
for doubling boundary: each of the 8 subregions we add the same number to all of the nodes and their children. WE just need to figure out the 8 numbers. 
for halving: same idea. if the subregion is getting 'deleted' then set its new tree index to UINT32_MAX (requires sorting? )

check if we need to do this after calculating forces (in calculate force method). for xthreaded, communicate the max dims between threads. 

double_and_sort
half_and_sort

maybe modify the sort method to do this?? will that work? i dont think so. we need two separate phases.  




{
    // Use IntelliSense to learn about possible attributes.
    // Hover to view descriptions of existing attributes.
    // For more information, visit: https://go.microsoft.com/fwlink/?linkid=830387
    "version": "2.0.0",
    "configurations": [
        
        {
            "name": "(gdb) Launch",
            "type": "cppdbg",
            "request": "launch",
            "preLaunchTask": "C/C++: g++.exe build active file",
            "program": "${fileDirname}/${fileBasenameNoExtension}.exe",
            "args": [],
            "stopAtEntry": false,
            "cwd": "${workspaceFolder}",
            "environment": [],
            "externalConsole": false,
            "MIMode": "gdb",
            "miDebuggerPath": "D:\\mingw-w64\\x86_64-8.1.0-posix-seh-rt_v6-rev0\\mingw64\\bin\\gdb.exe",
            "setupCommands": [
                {
                    "description": "Enable pretty-printing for gdb",
                    "text": "-enable-pretty-printing",
                    "ignoreFailures": true
                }
            ]
        },

    ]
}

*/
