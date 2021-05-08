#define __STDC_FORMAT_MACROS 1

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
#include <inttypes.h>


#define G (6.6743e-11) 
#define dt  0.1f
#define theta  1.0f
#define epsilon 0.1f
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

    uint32_t TAoS_index; // useless
} Particle;

typedef struct Node
{
    uint64_t tree_index; 
    
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
uint32_t N;

Particle* PAoS;
Node* TAoS;
Node* TAoS_swap_buffer; 

uint32_t TAoS_current_size;
uint32_t TAoS_max_size;  

float sys_boundary;
float tree_boundary; 


Index_Pair* sort_array;
uint32_t* map;


//float3* accels; // for debugging multithreading


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
        + "\ntotal mass: " +std::to_string(n.total_mass) +  "\n" + "PAoS index : " + 
            std::to_string(n.parent_TAoS_index) + "\nfirst_child_TAoS_index : " + std::to_string(n.first_child_TAoS_index) + "\n";  ;
}

int depth(int x) // N = 4 for quadtree, N = 8 for octree
{
    return  floor (log( (7) * x + 1) / log(8) );
}

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

bool is_internal_node(Node* node)
{
            // cant be empty node                     // must have children
    return(node->tree_index != UINT32_MAX &&  node->first_child_TAoS_index != UINT32_MAX);
}
bool is_blank_node(Node* node) // a blank is a leaf node thats not a particle. 
{
            // cant be empty node              // cant have children                         // can't have mass
    return (node->tree_index != UINT32_MAX &&  node->first_child_TAoS_index == UINT32_MAX && node->total_mass == 0.0f);
}
bool is_particle_node(Node* node)
{       
            // cant be empty node             // cant have children                          // mass must be nonzero
    return (node->tree_index != UINT32_MAX    && node->first_child_TAoS_index == UINT32_MAX && node->total_mass != 0.0f);
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

        n->tree_index = parent->tree_index*8 + i;
        n->first_child_TAoS_index = UINT32_MAX; // no children yet, these are blanks
        n->parent_TAoS_index = parent_TAoS_index;

        n->total_mass = 0.0f; // blank
        n->com = float3{0.0f, 0.0f, 0.0f};
    }
    //printf("inserted blanks as child of TAoS[%i]\n", parent_TAoS_index); // debug
}

void insert_particle_on_blank(Node* blank,  Node& particle_backup)
{
    if(!is_blank_node(blank))
    {
        printf("error in insert_particle_on_blank1 : node isn't blank!!!\n");
        return;
    }

    blank->com = particle_backup.com;
    blank->total_mass = particle_backup.total_mass;
    blank->first_child_TAoS_index = UINT32_MAX;
}

void insert_particle_on_blank(Node* blank,  uint32_t PAoS_index)
{
    if(!is_blank_node(blank))
    {
        printf("error in insert_particle_on_blank2 : node isn't blank!!!\n");
        return;
    }

    Particle* p =  &PAoS[PAoS_index];

    blank->com = p->pos;
    blank->total_mass = p->mass;
    blank->first_child_TAoS_index = UINT32_MAX;
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

    insert_blanks(TAoS_index);
    node->com = bkp.com;
    node->total_mass = bkp.total_mass;

    uint8_t sub_index = get_subregion_index(b, bkp.com);
    uint32_t target_TAoS_index = node->first_child_TAoS_index + sub_index; 

    if(!is_blank_node(&TAoS[target_TAoS_index]))
    {
        printf("TAoS[%i] ti %" PRIu64 "branch error at sub_index %i - internal on particle generated nonblank child at TAoS[%i]\n", TAoS_index, TAoS[TAoS_index].tree_index, sub_index, target_TAoS_index);
        printf("%s\n", nodeToString(TAoS[target_TAoS_index]).c_str());
    }


    // FIX THIS !!!
    insert_particle_on_blank(&TAoS[target_TAoS_index], bkp); // theres gotta be a better way than this. override the method to just take a Node& particle 
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

    //printf("merged new mass %6.4f\n", existing_particle->total_mass);
    //existing_particle->PAoS_index = PAoS_index; // can we just leave the PAoS index the same as the original? We cant have multiple values for this simultaneously.  
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


    // loop until we find a blank. branch tree if we hit a particle. 
    //while(!is_blank_node(current_node))  // while(is_internal) ? is that equivalent
    uint32_t depth = 0;
    uint32_t pinarow = 0;
    bool merged = false; 

    while(is_internal_node(current_node) || is_particle_node(current_node) && !merged)
    {   
        if(depth != 0)
        {
            parent_node = current_node;
        }
        depth++;
        
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
                //printf("\npinarow getting big : %i\n", pinarow);
                merge_particles(current_node, PAoS_index);
                merged = true;
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
    if(oob) { // safety check
        printf("PAoS[%i] was out of bounds and reached depth %i\n", PAoS_index, depth);
    }

    uint32_t parent_TAoS_index = (parent_node - TAoS); 
    if(is_blank_node(current_node))
    {
        insert_particle_on_blank(current_node, PAoS_index);
        return; 
    }
    /*else if (merged)
    {
        printf ("successfully merged particles:\n PAoS[%i] : (%6.4f, %6.4f, %6.4f)\n PAoS[%i] : (%6.4f, %6.4f, %6.4f)\n", current_node->PAoS_index, 
            PAoS[current_node->PAoS_index].pos.x, PAoS[current_node->PAoS_index].pos.y, PAoS[current_node->PAoS_index].pos.z, PAoS_index, PAoS[PAoS_index].pos.x, PAoS[PAoS_index].pos.y, PAoS[PAoS_index].pos.z);
        printf("new total mass: %6.4f\n", current_node->total_mass);
    }*/
    
}


void generate_TAoS(Node* TAoS, Particle* PAoS)
{
    TAoS_current_size = 0;

    // initialize root node
    Node* root = TAoS;
    root->tree_index = 0;
    root->parent_TAoS_index = UINT32_MAX; // no parent
    root->first_child_TAoS_index = 1; // we will put 8 blanks there
    root->com = float3{0.0f, 0.0f, 0.0f};
    root->total_mass = 0.0f; // not a particle

    
    TAoS_current_size++;

    // insert 8 blank children
    insert_blanks(0);

    // fill rest of array with 'empty nodes' 
    for(int ti = 9; ti < TAoS_max_size; ti++)
    {
        Node* n = &TAoS[ti];
        n->tree_index = UINT32_MAX;
        n->first_child_TAoS_index = UINT32_MAX;
        n->parent_TAoS_index = UINT32_MAX;
        
        n->total_mass = 0.0f; // new
        n->com = float3(0.0f, 0.0f, 0.0f);
    }

    for(int pi = 0; pi < N; pi++) {
        insert_particle_in_tree(TAoS, pi);
    }
}



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
    std::pair<Node*, float> start(TAoS, 2.0f*tree_boundary);
    stack.push(start); 

    // traverse
    while(!stack.empty())
    {
        std::pair<Node*, float> pair = stack.top();
        stack.pop();

        Node* current = pair.first;
        float width = pair.second;

        float3 r = current->com + (p->pos * -1.0f); 
        float r_norm = norm(r);

        ///printf("\t%i\tinternal: %i\t r_norm/width %6.2f\t total_mass %6.2f\t", current->tree_index, is_internal_node(current), r_norm / width, current->total_mass );

        if(is_internal_node(current) && (r_norm  / width) < theta ) // theta can be tuned
        {
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
            float3 a(0.0f, 0.0f, 0.0f);
            a = ((r / (r_norm + epsilon)) * (float) (G * current->total_mass)) / ((r_norm + epsilon ) * (r_norm + epsilon));
            accel =  a + accel;   
        }
    }

    p->vel = p->vel + accel * dt;
    p->pos = p->pos + p->vel*dt; 
    

    if(std::abs(p->pos.x) > sys_boundary)
        sys_boundary = std::abs(p->pos.x);
        
    if(std::abs(p->pos.y) > sys_boundary)
        sys_boundary = std::abs(p->pos.y);
        
    if(std::abs(p->pos.z) > sys_boundary)
        sys_boundary = std::abs(p->pos.z);
    
    tree_boundary = sys_boundary * boundary_fos; // only safe if we regenerate tree from scratch
    
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

#define show_percentages 1
int main (int argc, char** argv)
{
    uint8_t num_threads;
    uint32_t steps; 

    N = 1e3; 
    num_threads = 8;
    steps = 1e4; 


    if(argc < 4) {
        printf("\nUsage:\t%s <num_particles> <num_steps> <num_threads>\n", argv[0]);
        return -1; 
    }
    else {
        try {
            N = atoi(argv[1]);
            steps = atoi(argv[2]);
            num_threads = atoi(argv[3]);
        }
        catch(const char* msg)
        {
            printf("%s\n", msg);
            return -1;
        }
    }

    printf("\nSimulating %i particles over ", N);
    printf("%i timesteps of %4.2f s using %i threads...\n", steps, dt, num_threads);

    uint32_t *TMS, *TCS;
    TMS = &TAoS_max_size;
    TCS = &TAoS_current_size;


    //float3 sys_com = float3{0.0f, 0.0f, 0.0f};
    //float sys_mass = 0.0f;
    float max_dim = 0.0f;

    PAoS = (Particle*) malloc(N * sizeof(Particle)); 
    
    /*
    // SUN
    //PAoS[0].pos=float3(0.0f, 0.0f, 0.0f);
    float m_sun = 1.989e30;
    PAoS[0].pos=float3(-1.f, -1.f, -1.f); // slightly off center
    PAoS[0].vel = float3(0.0f, 0.0f, 0.0f);
    PAoS[0].mass = (m_sun); 
        

    // Small particle in orbit
    float r0 = 1.5e9;
    float v0 = (float) sqrt(G*m_sun / r0);
    PAoS[1].pos=float3(r0, 0.0f, 0.0f);   // r = 1.5*10^6 km
    PAoS[1].vel=float3(0.0f, v0, 0.0f); // v = 29.78 km/s              // Vorb = sqrt(G*Msun/R)
    PAoS[1].mass = 10000;                        // m = 1 kg

    */
    /* // EARTH
    PAoS[1].pos=float3(1.5e9, 0.0f, 0.0f);   // r = 1.5*10^6 km
    PAoS[1].vel=float3(0.0f, 2.978e4, 0.0f); // v = 29.78 km/s
    PAoS[1].mass = 5.972e24;                 // m = 5.972*10^24 kg
      */  

    
    
    srand(time(NULL));
    uint32_t range = 1000; 
    for(int i = 0; i < N; i++)
    {
        Particle* p = &PAoS[i];
        
        if(i >= 0) // here to optionally skip manually inserted particles
        {
            p->pos=float3( rand() % (range) - range/2.0f, rand() % range - range/2.0f, rand() % range - range/2.0f);
            //p->pos=float3(-500.0f+10.0f*i + rand()%10, -500.0f+10.0f*i +rand()%10, -500.0f+10.0f*i + rand()%10);
            p->vel=float3(0,0,0);
            p->mass = 1;
        }


        
        if(std::abs(p->pos.x) > max_dim)
            max_dim = std::abs(p->pos.x);
        
        if(std::abs(p->pos.y) > max_dim)
            max_dim = std::abs(p->pos.y);
        
        if(std::abs(p->pos.z) > max_dim)
            max_dim = std::abs(p->pos.z);

        //sys_com = sys_com + ((p->pos) * p->mass );
        //sys_mass += p->mass;
    }

    //sys_com = sys_com * (1.0f/sys_mass);
    sys_boundary = max_dim;
    tree_boundary = max_dim * boundary_fos;

    TAoS_max_size = fos * N;

    if(TAoS_max_size < 1000)
        TAoS_max_size = 1000; // necessary for small number of particles. ratio might get large
    
    TAoS = (Node*) malloc (TAoS_max_size * sizeof(Node));
    
    // we dont need this stuff unless we are sorting
    //TAoS_swap_buffer = (Node*) malloc (TAoS_max_size * sizeof(Node)); 
    //sort_array = (Index_Pair*) malloc(TAoS_max_size * sizeof(Index_Pair));
    //map = (uint32_t*) malloc(TAoS_max_size * sizeof(uint32_t));

    /*
    uint64_t update_TAoS_time = 0;
    uint64_t sort_TAoS_time = 0;
    uint64_t calculate_forces_time = 0;
    uint64_t gen_TAoS_time = 0;
    */
    uint64_t total = 0;


    auto t0 =  std::chrono::steady_clock::now();
    generate_TAoS(TAoS, PAoS);

    //gen_TAoS_time = (std::chrono::steady_clock::now() - t0).count();
    //uint32_t max_TAoS_size_reached = TAoS_current_size;
    


    //uint32_t gap = 100/dt; // every 100 seconds
    //float3* rSun2Earth = (float3*)malloc(steps / (gap/dt) * sizeof(float3) );
    auto t1 = std::chrono::steady_clock::now();
    for(int i = 0; i < steps; i++)
    {
        //auto t1 = std::chrono::steady_clock::now();
        generate_TAoS(TAoS, PAoS);
        //update_TAoS(TAoS, PAoS);
        //auto t2 = std::chrono::steady_clock::now();
        //sort_TAoS(TAoS, TAoS_swap_buffer, sort_array, map, TAoS_current_size, TAoS_max_size); // this is optional for cpu version. It just slows us down unecessarily! 
        //auto t3 = std::chrono::steady_clock::now();
        //calculate_forces(TAoS, PAoS, N);
        xthreaded_calculate_forces(num_threads, TAoS, PAoS, N);
        //auto t4 = std::chrono::steady_clock::now();

        /*
        if(TAoS_current_size > max_TAoS_size_reached)	
            max_TAoS_size_reached = TAoS_current_size;	

        update_TAoS_time += static_cast<uint64_t>((t2 - t1).count());
        sort_TAoS_time += static_cast<uint64_t>((t3 - t2).count());
        calculate_forces_time += static_cast<uint64_t>((t4 - t3).count());
        */
        
#if show_percentages
        if( i % (steps / 10) == 0) 
            printf("\t%i%%\n", (int) (10.f * (i / (steps / 10) ))   );
#endif

#if 0
        if(i % 10000 == 0)
        {
            float3 r(PAoS[1].pos.x - PAoS[0].pos.x, 
                    PAoS[1].pos.y - PAoS[0].pos.y,
                    PAoS[1].pos.z - PAoS[0].pos.z);

            

            //printf("\t(%3.1f,\t%3.1f,\t%3.1f)\t", r.x/1e6, r.y/1e6, r.z/1e6);
            printf("(%i,\t%i,\t%i)\t", (int) (r.x/1e6), (int)(r.y/1e6), (int)(r.z/1e6));

            printf("D=%i\t\t", (int) (sqrtf(r.x*r.x + r.y*r.y + r.z*r.z)/1e6));

            printf("@t=%is\n", (int)(i*dt));
        }
#endif
    }
    auto t2 = std::chrono::steady_clock::now();
    total = static_cast<uint64_t>((t2 - t1).count());

#if show_percentages
    printf("\t100%%\n");
#endif

    //update_TAoS_time /= steps;
    //sort_TAoS_time /= steps;
    //calculate_forces_time /= steps;
    
    //total = update_TAoS_time + sort_TAoS_time + calculate_forces_time;


    printf("\n\n");
    printf("PAoS size :\t\t%i\n", N);
    printf("final TAoS size :\t%i\n", TAoS_current_size);
    printf("ratio :\t\t\t%4.2f\n\n", ((TAoS_current_size * 1.0f) / N));

    //printf("max TAoS size :  %i\n", max_TAoS_size_reached);
    //printf("ratio     :  %4.2f\n", ((max_TAoS_size_reached * 1.0f) / N));

    /*
    printf("\nAvg Times:\n");
    printf("generate_TAoS_time : %lf ms\n", gen_TAoS_time/1E6);
    printf("update_TAoS_time : %lf ms\n", update_TAoS_time/1E6);
    printf("sort_TAoS_time : %lf ms\n", sort_TAoS_time/1E6);
    printf("calculate_forces_time : %lf ms\n", calculate_forces_time/1E6);
    */

    
    printf("total time :\t\t%6.3lf s\n\n", total/1E9); 
    printf("time per iteration :\t%6.3lf ms\n\n", total/1E6/steps);


    free(TAoS);
    free(TAoS_swap_buffer);
    free(PAoS);
    free(sort_array);
    free(map);

    return 1;

/*
   bounds b1{-1.0f, 1.0f, -1.0f, 1.0f, -1.0f, 1.0f};
   float3 pos = {-0.5f, 0.5f, 0.5f};

   printf("sub_index : %i\n\n", get_subregion_index(b1, pos));


   printf("index: %i ", (&TAoS[26] - TAoS) );
*/
}





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
*/


/*
profiling 
new solutions - toy example (similar to paper on slack)


bottleneck:
    as we introduce many threads, single threaded tree generation becomes a significant bottleneck

solution: 
    multithreaded tree generation 

        somewhat tricky bc traversal writes to nodes on the way down
            this would require atomic ops or mutex locking 


        we could make threads traversal in read only until they reach their position in the tree
        
        after building the tree, we then go back and add mass contributions.

        if we have 8 threads, this is simple. For each of the roots 8 subchildren, 
        one thread is responsible for that subtree. 

        Similarly simple to divide work for any # of threads thats a power of 2. 
        
        But they need to know which particles touched which nodes. when the threads traversed in read only mode, they needed to build a shared contributions list. 


    


*/



/*
PDF paper format:
YYYY-where published-title
*/