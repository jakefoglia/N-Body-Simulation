using UnityEngine;
using System.Collections;
using System.Collections.Generic;
using System;

using UnityEditor;

/*
 * TODOs
 * 
 *      
 * 
 * */

 static class Constants {
    public const float G = 10000;//6.6743e-11f;
    public const float dt = 0.25f;//0.1f;
    public const float theta = 0.75f; // higher theta means more accurate, lower means faster
    public const float epsilon = 0.5f; 
    public const int fos = 20;
    public const float boundary_fos = 1.25f;
    public const uint N = 1000u;
}




public class ParticleManager : MonoBehaviour
{
    // global arrays and vars
    bool simulationOn;
    //Particle[] PAoS;
    Node[] TAoS;
    Node[] TAoS_swap_buffer;
    GameObject[] particle_objects;

    uint TAoS_current_size;
    uint TAoS_max_size;

    float sys_boundary;
    float tree_boundary;

    Index_Pair[] sort_array;
    uint[] map;

    [SerializeField] GameObject particle_prefab;

    Gradient gradient;
    GradientColorKey[] colorKey;
    GradientAlphaKey[] alphaKey;

    Color color_purple;
    Color color_orange;

    // Start is called before the first frame update
    void Start()
    {
        //num_threads = 8;
        //steps = 10;// (int) 1e4;


        //Vector3 sys_com = Vector3{0.0f, 0.0f, 0.0f};
        //float sys_mass = 0.0f;

        //PAoS = (Particle*)malloc(N * sizeof(Particle));

        // populate and init particles
        init(ref TAoS, ref particle_objects, InitializationMode.RandomUniformSphere, ref sys_boundary, ref tree_boundary, Constants.N, 1000);


        //srand(time(NULL));

        //sys_com = sys_com * (1.0f/sys_mass);

        var cam = UnityEngine.Object.FindObjectOfType<Camera>();
        cam.clearFlags = CameraClearFlags.SolidColor;
        cam.backgroundColor = Color.black;

        gradient = new Gradient();

        // Populate the color keys at the relative time 0 and 1 (0 and 100%)
        color_purple  = new Color(105f / 255f, 0f / 255f, 194f / 255f);
        color_orange = new Color(203f / 255f, 8f / 255f, 8f / 255f);

        var cold = color_purple;
        var hot = color_orange;

        colorKey = new GradientColorKey[2];
        colorKey[0].color = cold;
        colorKey[0].time = 0.0f;
        colorKey[1].color = hot;
        colorKey[1].time = 1.0f;

        // Populate the alpha  keys at relative time 0 and 1  (0 and 100%)
        alphaKey = new GradientAlphaKey[2];
        alphaKey[0].alpha = 1.0f;
        alphaKey[0].time = 0.0f;
        alphaKey[1].alpha = 0.0f;
        alphaKey[1].time = 1.0f;

        gradient.SetKeys(colorKey, alphaKey);
    }

    // Update is called once per frame
    void Update()
    {
        generate_TAoS(ref TAoS, ref particle_objects, Constants.N); 
   
        calculate_forces(ref TAoS, ref particle_objects, Constants.N, Time.deltaTime);
        //xthreaded_calculate_forces(ref TAoS, ref PAoS, Constants.N, Time.deltaTime);

        //update_particle_objects(ref PAoS, ref particle_objects, Constants.N); 
    }

    enum InitializationMode
    {
        RandomUniformSphere = 0,
        SolarSystem = 1
    }

    void init(ref Node[] TAoS, ref GameObject[] particle_objects, InitializationMode initMode, ref float sys_boundary, ref float tree_boundary, uint N = 0, float range = 1000)
    {
        float max_dim = 0.0f;

        switch (initMode)
        {
            case InitializationMode.RandomUniformSphere:
                {
                    TAoS_max_size = Constants.fos * N; 
                    
                    if (TAoS_max_size < 1000)
                        TAoS_max_size = 1000; // necessary for small number of particles. ratio might get large

                    TAoS = new Node[TAoS_max_size];

                    particle_objects = new GameObject[N];

                    for (uint i = 0; i < N; i++)
                    {
                        particle_objects[i] = Instantiate(particle_prefab) as GameObject;
                    }

                    for (int i = 0; i < N; i++)
                    {
                        ref GameObject p = ref particle_objects[i];

                        var u = UnityEngine.Random.Range(0f, 1f);
                        var v = UnityEngine.Random.Range(0f, 1f);
                        float r = UnityEngine.Random.Range(0f, 1f); 

                        float angle1 = 2 * Mathf.PI * u;
                        float angle2 = 2 * Mathf.PI * v;
                        //float angle2 = Mathf.Pow(Mathf.Cos(2 * v - 1), -1);

                        float X = (r * range * Mathf.Cos(angle1) * Mathf.Cos(angle2));
                        float Y = (r * range * Mathf.Sin(angle1) * Mathf.Cos(angle2));
                        float Z = (r * range * Mathf.Sin(angle2));
                        
                        p.transform.localPosition = new Vector3(X, Y, Z);
                        //p.GetComponent<Rigidbody>().velocity = Vector3.zero;

                        p.GetComponent<Rigidbody>().velocity = Vector3.zero;
                        p.GetComponent<Rigidbody>().mass = 1;                     

                        if (Mathf.Abs(p.transform.localPosition.x) > max_dim)
                            max_dim = Mathf.Abs(p.transform.localPosition.x);

                        if (Mathf.Abs(p.transform.localPosition.y) > max_dim)
                            max_dim = Mathf.Abs(p.transform.localPosition.y);

                        if (Mathf.Abs(p.transform.localPosition.z) > max_dim)
                            max_dim = Mathf.Abs(p.transform.localPosition.z);

                        //sys_com = sys_com + ((p.transform.localPosition) * p.GetComponent<Rigidbody>().mass );
                        //sys_mass += p.GetComponent<Rigidbody>().mass;
                    }

                    sys_boundary = max_dim;
                    tree_boundary = max_dim * Constants.boundary_fos;

                break;
            }
            case InitializationMode.SolarSystem:
                {
                    //PAoS = new Particle[9];
                    /*
                    // SUN
                    //particle_objects[0].pos=Vector3.zero;
                    float m_sun = 1.989e30;
                    particle_objects[0].pos=new Vector3(-1.f, -1.f, -1.f); // slightly off center
                    particle_objects[0].vel = Vector3.zero;
                    particle_objects[0].mass = (m_sun); 

                    // Small particle in orbit
                    float r0 = 1.5e9;
                    float v0 = (float) sqrt(G*m_sun / r0);
                    particle_objects[1].pos=new Vector3(r0, 0.0f, 0.0f);   // r = 1.5*10^6 km
                    particle_objects[1].vel=new Vector3(0.0f, v0, 0.0f); // v = 29.78 km/s              // Vorb = sqrt(G*Msun/R)
                    particle_objects[1].mass = 10000;                        // m = 1 kg

                    */
                    /* // EARTH
                    particle_objects[1].pos=new Vector3(1.5e9, 0.0f, 0.0f);   // r = 1.5*10^6 km
                    particle_objects[1].vel=new Vector3(0.0f, 2.978e4, 0.0f); // v = 29.78 km/s
                    particle_objects[1].mass = 5.972e24;                 // m = 5.972*10^24 kg
                      */

                    for (int i = 0; i < N; i++)
                    {
                        ref GameObject p = ref particle_objects[i];

                        if (Mathf.Abs(p.transform.localPosition.x) > max_dim)
                            max_dim = Mathf.Abs(p.transform.localPosition.x);

                        if (Mathf.Abs(p.transform.localPosition.y) > max_dim)
                            max_dim = Mathf.Abs(p.transform.localPosition.y);

                        if (Mathf.Abs(p.transform.localPosition.z) > max_dim)
                            max_dim = Mathf.Abs(p.transform.localPosition.z);

                        //sys_com = sys_com + ((p.transform.localPosition) * p.GetComponent<Rigidbody>().mass );
                        //sys_mass += p.GetComponent<Rigidbody>().mass;
                    }

                    sys_boundary = max_dim;
                    tree_boundary = max_dim * Constants.boundary_fos;

                    break;
            }
        }
    }

    Color colorMap(float v) // v should be positive (magnitude of either vel or acceleration)
    {
        const float vFloor = 0f;
        //const float vCeil = 150f; // good for using vel
        const float vCeil = 500f; // good for using accel   todo calculate a good range based on mass and grav 

        v = Mathf.Clamp((v - vFloor) / (vCeil - vFloor), 0f, 1f);
        return gradient.Evaluate(v);
    }
    
    void update_particle_objects(ref ParticleCollisionEvent[] PAoS, ref GameObject[] particle_objects, uint N)
    {
        //for(uint i = 0; i < n; i++) { 
        //    particle_objects[i].transform.localPosition = particle_objects[i].pos;
            /*
            var renderer = particle_objects[i].gameObject.GetComponent<Renderer>();
            renderer.material.SetColor("_Color", colorMap(particle_objects[i].vel.magnitude)); // were gonna use accel 
            */
        //}
    }

    


    
    string nodeToString(Node n)
    {
        string type = "";
        if (is_empty_node(n))
            type += "empty ";
        if (is_blank_node(n))
            type += "blank ";
        if (is_internal_node(n))
            type += "internal ";
        if (is_particle_node(n))
            type += "particle ";

        return ("tree index : ") + (n.tree_index) + "  " + type + "\ncom: (" + (n.com.x) + ", " + (n.com.y) + ", " + (n.com.z) + ")"
            + "\ntotal mass: " + (n.total_mass) + "\n" + "PAoS index : " +
                (n.parent_TAoS_index) + "\nfirst_child_TAoS_index : " + (n.first_child_TAoS_index) + "\n"; ;
    }

    int depth(int x) // N = 4 for quadtree, N = 8 for octree
    {
        return Mathf.FloorToInt(Mathf.Log(7 * x + 1, 8));
    }

    byte get_subregion_index(bounds parent_region, Vector3 pos) // returns 0 through 7, specifying the subregion. x is 3rd bit, y is 2nd bit, z is 1st bit
    {
        Vector3 center;
        center.x = (parent_region.x_lb + parent_region.x_ub) / 2.0f;
        center.y = (parent_region.y_lb + parent_region.y_ub) / 2.0f;
        center.z = (parent_region.z_lb + parent_region.z_ub) / 2.0f;

        int index = 0;
        if (pos.x >= center.x)
            index = index | (1 << 2);
        if (pos.y >= center.y)
            index = index | (1 << 1);
        if (pos.z >= center.z)
            index = index | (1 << 0);

        return (byte)   index;
    }
    bounds get_subregion(bounds parent_region, byte index)
    {
        Vector3 center;
        center.x = (parent_region.x_lb + parent_region.x_ub) / 2.0f;
        center.y = (parent_region.y_lb + parent_region.y_ub) / 2.0f;
        center.z = (parent_region.z_lb + parent_region.z_ub) / 2.0f;
        bounds region;
        if ((index & (1 << 2)) != 0) // x upper half
        {
            region.x_lb = center.x;
            region.x_ub = parent_region.x_ub;
        }
        else
        {
            region.x_lb = parent_region.x_lb;
            region.x_ub = center.x;
        }
        if ((index & (1 << 1)) != 0) // y upper half
        {
            region.y_lb = center.y;
            region.y_ub = parent_region.y_ub;
        }
        else
        {
            region.y_lb = parent_region.y_lb;
            region.y_ub = center.y;
        }
        if ((index & (1 << 2)) != 0) // z upper half
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
    bounds get_subregion(bounds parent_region, Vector3 pos)
    {
        //byte index = get_subregion_index(parent_region, pos);
        //return get_subregion(parent_region, index);

        Vector3 center;
        center.x = (parent_region.x_lb + parent_region.x_ub) / 2.0f;
        center.y = (parent_region.y_lb + parent_region.y_ub) / 2.0f;
        center.z = (parent_region.z_lb + parent_region.z_ub) / 2.0f;

        bounds region;
        if (pos.x >= center.x)
        {
            region.x_lb = center.x;
            region.x_ub = parent_region.x_ub;
        }
        else
        {
            region.x_lb = parent_region.x_lb;
            region.x_ub = center.x;
        }
        if (pos.y >= center.y)
        {
            region.y_lb = center.y;
            region.y_ub = parent_region.y_ub;
        }
        else
        {
            region.y_lb = parent_region.y_lb;
            region.y_ub = center.y;
        }
        if (pos.z >= center.z)
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

    bool is_child_tree_index(uint pti, uint cti)
    {
        return ((cti - 1) / 8 == pti);
    }

    bool is_internal_node(Node node)
    {
        // cant be empty node                     // must have children
        return (node.tree_index != UInt32.MaxValue && node.first_child_TAoS_index != UInt32.MaxValue);
    }
    bool is_blank_node(Node node) // a blank is a leaf node thats not a particle. 
    {
        // cant be empty node              // cant have children                         // can't have mass
        return (node.tree_index != UInt32.MaxValue && node.first_child_TAoS_index == UInt32.MaxValue && node.total_mass == 0.0f);
    }
    bool is_particle_node(Node node)
    {
        // cant be empty node             // cant have children                          // mass must be nonzero
        return (node.tree_index != UInt32.MaxValue && node.first_child_TAoS_index == UInt32.MaxValue && node.total_mass != 0.0f);
    }
    bool is_empty_node(Node node)
    {
        // all nonempty nodes must have tree index
        return (node.tree_index == UInt32.MaxValue); // index out of bounds means we treat as empty
    }


    // 3 types of nodes: internal, particle, blank. 
    // root starts as blank: it has no children but its not a particle
    // then after trying to add one particle, we add the particle and 7 other blank node children which will be turned into particles then internal nodes when necessary. 

    void insert_blanks(uint parent_TAoS_index)
    {
        ref Node parent = ref  TAoS[parent_TAoS_index];
        parent.first_child_TAoS_index = TAoS_current_size;

        // insert 8 blank children
        for (int i = 1; i <= 8; i++)
        {
            ref Node n = ref TAoS[TAoS_current_size++];

            ulong v = (parent.tree_index * 8);
            n.tree_index = v + (ulong) i;
            n.first_child_TAoS_index = UInt32.MaxValue; // no children yet, these are blanks
            n.parent_TAoS_index = parent_TAoS_index;

            n.total_mass = 0.0f; // blank
            n.com = new Vector3( 0.0f, 0.0f, 0.0f);
        }
        //Debug.Log("inserted blanks as child of TAoS[%i]\n", parent_TAoS_index); // debug
    }

    void insert_particle_on_blank(ref Node blank, ref Node particle_backup)
    {
        if (!is_blank_node(blank))
        {
            Debug.Log("error in insert_particle_on_blank1 : node isn't blank!!!\n");
            return;
        }

        blank.com = particle_backup.com;
        blank.total_mass = particle_backup.total_mass;
        blank.first_child_TAoS_index = UInt32.MaxValue;
    }

    void insert_particle_on_blank(ref Node blank, uint PAoS_index)
    {
        if (!is_blank_node(blank))
        {
            Debug.Log("error in insert_particle_on_blank2 : node isn't blank!!!\n");
            return;
        }

        ref GameObject p = ref particle_objects[PAoS_index];

        blank.com = p.transform.localPosition;
        blank.total_mass = p.GetComponent<Rigidbody>().mass;// p.GetComponent<Rigidbody>().mass;
        blank.first_child_TAoS_index = UInt32.MaxValue;
    }

    void insert_internal_on_particle(uint TAoS_index, bounds b)
    {
        ref Node node = ref TAoS[TAoS_index];
        if (!is_particle_node(node))
        {
            Debug.Log("error in insert_inernal_on_particle : node isn't particle!!!\n");
            return;
        }

        Node bkp = node.clone(); // create a copy

        insert_blanks(TAoS_index);
        node.com = bkp.com;
        node.total_mass = bkp.total_mass;

        byte sub_index = get_subregion_index(b, bkp.com);
        uint target_TAoS_index = node.first_child_TAoS_index + sub_index;

        if (!is_blank_node(TAoS[target_TAoS_index]))
        {
            Debug.Log(String.Format("TAoS[%i] ti %", "branch error at sub_index %i - internal on particle generated nonblank child at TAoS[%i]\n", TAoS_index, TAoS[TAoS_index].tree_index, sub_index, target_TAoS_index));
            Debug.Log(String.Format("%s\n", nodeToString(TAoS[target_TAoS_index])));
        }


        // FIX THIS !!!
        insert_particle_on_blank(ref TAoS[target_TAoS_index], ref bkp); // theres gotta be a better way than this. override the method to just take a Node& particle 
    }

    void merge_particles(Node existing_particle, uint PAoS_index)
    {
        if (!is_particle_node(existing_particle))
        {
            Debug.Log("error in merge_particles : existing node isn't a particle!!!\n");
            return;
        }

        ref GameObject p = ref particle_objects[PAoS_index];

        existing_particle.com = ((existing_particle.com * existing_particle.total_mass) + p.transform.localPosition * p.GetComponent<Rigidbody>().mass) / (existing_particle.total_mass + p.GetComponent<Rigidbody>().mass);
        existing_particle.total_mass = (existing_particle.total_mass + p.GetComponent<Rigidbody>().mass);

        //Debug.Log("merged new mass %6.4f\n", existing_particle.total_mass);
        //existing_particle.PAoS_index = PAoS_index; // can we just leave the PAoS index the same as the original? We cant have multiple values for this simultaneously.  
    }

    void insert_particle_in_tree(ref Node[] TAoS, uint PAoS_index) // a node should have 0 children, or 8 children (must be consecutive)
    {
        //Debug.Log("\ninserting particle_objects[%i] TCS: %i  TMS: %i  POS(%6.4f, %6.4f, %6.4f)\n", PAoS_index, TAoS_current_size, TAoS_max_size, particle_objects[PAoS_index].pos.x, particle_objects[PAoS_index].pos.y, particle_objects[PAoS_index].pos.z);
        ref GameObject p = ref particle_objects[PAoS_index];
        bounds current_bounds = new bounds(-tree_boundary, tree_boundary, -tree_boundary, tree_boundary, -tree_boundary, tree_boundary); // start with entire tree boundary
        uint current_node_index = 0;
        uint parent_node_index = 0;


        bool oob = false;
        if (particle_objects[PAoS_index].transform.localPosition.x < current_bounds.x_lb || particle_objects[PAoS_index].transform.localPosition.x > current_bounds.x_ub ||
        particle_objects[PAoS_index].transform.localPosition.y < current_bounds.y_lb || particle_objects[PAoS_index].transform.localPosition.y > current_bounds.y_ub ||
        particle_objects[PAoS_index].transform.localPosition.z < current_bounds.z_lb || particle_objects[PAoS_index].transform.localPosition.z > current_bounds.z_ub)
        {
            Debug.Log("PARTICLE IS OUT OF BOUNDS!!!!!!!!!\n");
            Debug.Log(String.Format("(%6.4f, %6.4f, %6.4f)\t(%6.4f, %6.4f, %6.4f) - (%6.4f, %6.4f, %6.4f)\n", particle_objects[PAoS_index].transform.localPosition.x, particle_objects[PAoS_index].transform.localPosition.y, particle_objects[PAoS_index].transform.localPosition.z, current_bounds.x_lb, current_bounds.y_lb, current_bounds.z_lb,
                            current_bounds.x_ub, current_bounds.y_ub, current_bounds.z_ub));
            oob = true;
        }


        // loop until we find a blank. branch tree if we hit a particle. 
        //while(!is_blank_node(current_node))  // while(is_internal) ? is that equivalent
        uint depth = 0;
        uint pinarow = 0;
        bool merged = false;

        while (is_internal_node(TAoS[current_node_index]) || is_particle_node(TAoS[current_node_index]) && !merged)
        {
            if (depth != 0)
            {
                parent_node_index = current_node_index;
            }
            depth++;

            byte sub_index = get_subregion_index(current_bounds, p.transform.localPosition);
            current_bounds = get_subregion(current_bounds, sub_index);

            if (TAoS[current_node_index].first_child_TAoS_index == 0)
                Debug.Log("uhh ohh 1 in insert_node()\n\n");
            if ((!is_particle_node(TAoS[current_node_index]) && TAoS[current_node_index].first_child_TAoS_index == UInt32.MaxValue)) // debug
                Debug.Log("uhh ohh 2 in insert_node()\n\n");

            if (is_particle_node(TAoS[current_node_index]))
            {
                pinarow++;
                if (pinarow == 30) // problem is that two particle have very similar positions. just combine them into a single particle node on the tree. 
                {
                    //Debug.Log("\npinarow getting big : %i\n", pinarow);
                    merge_particles(TAoS[current_node_index], PAoS_index);
                    merged = true;
                }
                else
                    insert_internal_on_particle(current_node_index, current_bounds);
            }
            else // only advance if it was already and internal node. 
            {
                // add mass contribution
                TAoS[current_node_index].com = ((TAoS[current_node_index].com * TAoS[current_node_index].total_mass) + (p.transform.localPosition * p.GetComponent<Rigidbody>().mass)) * (1.0f / (TAoS[current_node_index].total_mass + p.GetComponent<Rigidbody>().mass));
                TAoS[current_node_index].total_mass += p.GetComponent<Rigidbody>().mass;

                // this is why we must insert 8 children at once, or none 
                current_node_index = TAoS[current_node_index].first_child_TAoS_index + sub_index;
            }
        }
        if (oob)
        { // safety check
            Debug.Log(String.Format("particle_objects[%i] was out of bounds and reached depth %i\n", PAoS_index, depth));
        }

        if (is_blank_node(TAoS[current_node_index]))
        {
            insert_particle_on_blank(ref TAoS[current_node_index], PAoS_index);
            return;
        }
        /*else if (merged)
        {
            printf ("successfully merged particles:\n particle_objects[%i] : (%6.4f, %6.4f, %6.4f)\n particle_objects[%i] : (%6.4f, %6.4f, %6.4f)\n", current_node.PAoS_index, 
                particle_objects[current_node.PAoS_index].pos.x, particle_objects[current_node.PAoS_index].pos.y, particle_objects[current_node.PAoS_index].pos.z, PAoS_index, particle_objects[PAoS_index].pos.x, particle_objects[PAoS_index].pos.y, particle_objects[PAoS_index].pos.z);
            Debug.Log("new total mass: %6.4f\n", current_node.total_mass);
        }*/

    }


    void generate_TAoS(ref Node[] TAoS, ref GameObject[] particle_objects, uint N)
    {
        TAoS_current_size = 0;

        // initialize root node
        ref Node root = ref TAoS[0];
        root.tree_index = 0;
        root.parent_TAoS_index = UInt32.MaxValue; // no parent
        root.first_child_TAoS_index = 1; // we will put 8 blanks there
        root.com = new Vector3( 0.0f, 0.0f, 0.0f);
        root.total_mass = 0.0f; // not a particle


        TAoS_current_size++;

        // insert 8 blank children
        insert_blanks(0);


        // fill rest of array with 'empty nodes' 
        for (int ti = 9; ti < TAoS_max_size; ti++)
        {
            ref Node n = ref TAoS[ti];
            n.tree_index = UInt32.MaxValue;
            n.first_child_TAoS_index = UInt32.MaxValue;
            n.parent_TAoS_index = UInt32.MaxValue;

            n.total_mass = 0.0f; // new
            n.com = Vector3.zero;
        }

        for (uint pi = 0; pi < N; pi++)
        {
            insert_particle_in_tree(ref TAoS, pi);
        }
    }



    void sort_TAoS(ref Node[] TAoS, ref Node[] TAoS_swap_buffer, ref Index_Pair[] n2o, ref uint[] m, uint N, uint max)
    {
        for (uint i = 0; i < N; i++)
        {
            n2o[i].index1 = TAoS[i].tree_index;
            n2o[i].index2 = i; // TAoS_index_old
        }

        // insertion sort
        Index_Pair key;
        for (uint i = 1, j = 0; i < N; i++)
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
        for (uint i = 0; i < N; i++)
        {
            m[n2o[i].index2] = i; // map [TAoS_index_old] = TAoS_index_new
        }

        // root is special case ( has no map for parent_TAoS_index)
        TAoS_swap_buffer[0] = TAoS[0]; // todo unity this good? 
        TAoS_swap_buffer[0].first_child_TAoS_index = m[TAoS[0].first_child_TAoS_index];

        // remap remaining non-empty nodes
        for (int i = 1; i < N; i++)
        {
            TAoS_swap_buffer[m[i]] = TAoS[i];
            if (TAoS[i].first_child_TAoS_index < UInt32.MaxValue)
                TAoS_swap_buffer[m[i]].first_child_TAoS_index = m[TAoS[i].first_child_TAoS_index]; // if this seg faults, we are trying to map based on an invalid index

            if (TAoS[i].parent_TAoS_index < UInt32.MaxValue)
                TAoS_swap_buffer[m[i]].parent_TAoS_index = m[TAoS[i].parent_TAoS_index]; // if this seg faults, we are trying to map based on an invalid index
        }
        for (uint i = N; i < TAoS_max_size; i++)
        {
            TAoS_swap_buffer[i] = TAoS[i];
        }

        // swap the buffers // todo unity this good ?
        Node[] temp = TAoS;
        TAoS = TAoS_swap_buffer;
        TAoS_swap_buffer = temp;
    }

    float norm(Vector3 vec)
    {
        return Mathf.Sqrt(vec.x * vec.x + vec.y * vec.y + vec.z * vec.z);
    }


    void dfs_traverse(ref Node[] TAoS, uint particle_index, float frame_time)
    {
        ref GameObject p = ref particle_objects[particle_index];

        ///Debug.Log("\n--------------------------------------------------------------\n");
        ///Debug.Log("%i  (%6.4f, %6.4f, %6.4f)\n",p - PAoS, p.pos.x, p.pos.y, p.pos.z);

        Vector3 accel = Vector3.zero;
        float magnitudeSum = 0f;

        // init stack with root
        Stack<KeyValuePair<Node, float>> stack = new Stack<KeyValuePair<Node, float>>();
        KeyValuePair<Node, float> start = new KeyValuePair<Node, float> (TAoS[0], 2.0f * tree_boundary);
        stack.Push(start);

        // traverse
        while (stack.Count!=0)
        {
            KeyValuePair<Node, float> pair = stack.Pop();

            Node current = pair.Key;
            float width = pair.Value;

            Vector3 r = current.com + (p.transform.localPosition * -1.0f);
            float r_norm = norm(r);

            ///Debug.Log("\t%i\tinternal: %i\t r_norm/width %6.2f\t total_mass %6.2f\t", current.tree_index, is_internal_node(current), r_norm / width, current.total_mass );

            if (is_internal_node(current) && (r_norm / width) < Constants.theta) // theta can be tuned
            {
                // push valid children to stack
                ref Node child = ref TAoS[current.first_child_TAoS_index];
                for (int i = 0; i < 8; i++)
                {
                    child = ref TAoS[current.first_child_TAoS_index + i];
                    if (!is_empty_node(child) && !is_blank_node(child))
                        stack.Push(new KeyValuePair<Node, float>(child, width / 2.0f));
                }
            }
            else
            {
                Vector3 a;
                a = r;
                a *= (Constants.G * current.total_mass);
                a /= (r_norm + Constants.epsilon);
                a /= (r_norm + Constants.epsilon);
                a /= (r_norm + Constants.epsilon);
                accel += a;
                magnitudeSum += a.magnitude;
            }
        }


        //p.transform.localPosition = p.transform.localPosition + p.GetComponent<Rigidbody>().velocity * frame_time * Constants.dt; // handles by velocity
        var newVel = p.GetComponent<Rigidbody>().velocity + accel * frame_time * Constants.dt;
        p.GetComponent<Rigidbody>().velocity = newVel; // todo add force 
        


        var renderer = p.GetComponent<Renderer>(); // update color based on accel
        var color = colorMap(magnitudeSum);
        renderer.material.SetColor("_Color", color);
        var mat = p.GetComponent<Renderer>().material;
        mat.SetColor("_EmissionColor", color);



        if (Mathf.Abs(p.transform.localPosition.x) > sys_boundary)
            sys_boundary = Mathf.Abs(p.transform.localPosition.x);

        if (Mathf.Abs(p.transform.localPosition.y) > sys_boundary)
            sys_boundary = Mathf.Abs(p.transform.localPosition.y);

        if (Mathf.Abs(p.transform.localPosition.z) > sys_boundary)
            sys_boundary = Mathf.Abs(p.transform.localPosition.z);

        tree_boundary = sys_boundary * Constants.boundary_fos; // only safe if we regenerate tree from scratch

    }
    void calculate_forces(ref Node[] TAoS, ref GameObject[] particle_objects, uint N, float frame_time)
    {
        for (uint i = 0; i < N; i++)
        {
            dfs_traverse(ref TAoS, i, frame_time); // traverse starting at the root node
        }
    }
    /*
    void* xthreaded_traversal(void* arg)
    {
        xthread_arg_struct* arg_struct = (xthread_arg_struct*)arg;

        Node TAoS = arg_struct.TAoS;
        Particle PAoS_start_address = arg_struct.PAoS_start_address;
        uint N = arg_struct.n;

        calculate_forces(TAoS, PAoS_start_address, n);

        pthread_exit(0);
        return nullptr;
    }
    */
    /*

    void xthreaded_calculate_forces(byte num_threads, Node TAoS, Particle PAoS, uint N)
    {
        xthread_arg_struct arg_structs[num_threads];
        if (num_threads < 2)
        {
            calculate_forces(TAoS, PAoS, n);
        }
        else
        {
            pthread_t tids[num_threads - 1];
            Particle running_PAoS_ptr = PAoS;
            uint ppt = n / (num_threads - 1); // particle per thread
            for (int i = 0; i < num_threads - 1; i++)
            {
                arg_structs[i].TAoS = TAoS;
                arg_structs[i].PAoS_start_address = running_PAoS_ptr;
                arg_structs[i].n = ppt;

                pthread_attr_t attr;
                pthread_attr_init(&attr);

                pthread_create(&tids[i], &attr, xthreaded_traversal, (void*)&arg_structs[i]);
                running_PAoS_ptr = &running_PAoS_ptr[ppt];
            }
            calculate_forces(TAoS, running_PAoS_ptr, n - ppt * (num_threads - 1)); // traverse the remainder in the meantime;

            for (int i = 0; i < num_threads - 1; i++)
            {
                pthread_join(tids[i], NULL);
            }
        }
        return;
    }
    */
}


struct bounds
{
    public float x_lb, x_ub, y_lb, y_ub, z_lb, z_ub;
    public bounds(float x_lb, float x_ub, float y_lb, float y_ub, float z_lb, float z_ub) {
            this.x_lb = x_lb;
            this.x_ub = x_ub;
            this.y_lb = y_lb;
            this.y_ub = y_ub;
            this.z_lb = z_lb;
            this.z_ub = z_ub;
     }

    //inline bounds operator*(float s) const { return bounds(x_lb* s, x_ub* s, y_lb* s, y_ub* s, z_lb* s, z_ub* s);
}

/*
struct Particle
{
   public Vector3 pos;
   public Vector3 vel;
   public float mass;

   public uint TAoS_index; // useless
}

*/
struct Node
{
    public ulong tree_index;

    public uint parent_TAoS_index;
    public uint first_child_TAoS_index;

    public Vector3 com;
    public float total_mass;

    public Node clone()
    {
        Node n = new Node();
        n.tree_index = this.tree_index;
        n.parent_TAoS_index = this.parent_TAoS_index;
        n.first_child_TAoS_index = this.first_child_TAoS_index;
        n.com = this.com;
        n.total_mass = this.total_mass;

        return n;
    }
}

struct Index_Pair
{
    public ulong index1, index2;
}

/*
struct xthread_arg_struct
{
    Node TAoS;
    Particle PAoS_start_address;
    uint N;
}
*/

//Vector3* accels; // for debugging multithreading

/*
bool is_internal_node(Node node);
bool is_particle_node(Node node);
bool is_blank_node(Node node);
bool is_empty_node(Node node);
*/


/*
Idea for handling out of bounds without regenerating tree (dynamic adjustment)

when particle is oob, rather than branching infinitely when another oob particle competes with it, we just maintain a linked list for those particles.
we cam't really do the hashing strategy bc then we will unjustly occupy  positions for valid particles (in bounds).

oob linked list! or array ? 
maybe we double tree boundary if the array gets filled, in addition to sys boundary doubling or halving 


remap of tree indices is easier than you think. 
for doubling boundary: each of the 8 subregions we add the same number to all of the nodes and their children. WE just need to figure out the 8 numbers. 
for halving: same idea. if the subregion is getting 'deleted' then set its new tree index to UInt32.MaxValue (requires sorting? )

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

