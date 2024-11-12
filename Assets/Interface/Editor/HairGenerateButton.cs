using UnityEngine;
using UnityEditor;
using Unity.VisualScripting;

[CustomEditor(typeof(HairComponent))]
public class HairGenerateButton : Editor
{
    public override void OnInspectorGUI()
    {
        base.OnInspectorGUI();

        HairComponent component = (HairComponent)target;
        if (GUILayout.Button("Generate Hair Mesh"))
        {
            component.WhenGenerationButtonClicked();
        }
    }
}