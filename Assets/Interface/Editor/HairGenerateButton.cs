using UnityEngine;
using UnityEditor;
using Unity.VisualScripting;

[CustomEditor(typeof(HairComponent))]
public class HairGenerateButton : Editor
{

    private SerializedProperty _strandHairs;

   private void OnEnable()
    {
        _strandHairs = serializedObject.FindProperty("strandHairs");
    }
    public override void OnInspectorGUI()
    {
        base.OnInspectorGUI();

        HairComponent component = (HairComponent)target;
        if (GUILayout.Button("Generate Hair Mesh"))
        {
            if (component.GetNumOfHairPram() > 0)
            {
                foreach(var Param in _strandHairs)
                {
                    //addfur
                    Debug.Log("Hair Mesh" + /* index*/ "Created!");
                }
            }
            else
            {
                Debug.LogError("There's no StrandHairParameter Instance!");
            }
        }
    }
}